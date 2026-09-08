using Test
using LinearAlgebra
using MatrixPade

isdefined(@__MODULE__, :_torder_coeff) || include("matrix_pade_testutils.jl")
isdefined(@__MODULE__, :_texp_iab_series) || include("exp_series_testutils.jl")

# The README's worked example: f(x) = exp(i A x + i B x^2) with
#
#   A = [1 i; -i -1],  B = [0 1; 1 0].
#
# A and B anticommute, so the eigenbasis of M(x) = i(A x + B x^2) depends on x
# and the series cannot be diagonalised once and approximated eigenvalue by
# eigenvalue. The series coefficients below are the ones printed by
# Mathematica's Series[MatrixExp[I A x + I B x^2], {x, 0, 12}].
@testset "README example: exp(i A x + i B x^2)" begin
    T = Complex{Rational{BigInt}}
    A = _treadme_A(T)
    B = _treadme_B(T)
    coeffs = _texp_iab_series(A, B, 12)

    @testset "series matches the published expansion" begin
        e11 = T[1, im, -1, -im//3, -1//3, -2im//15, 7//45, 2im//63, 8//315,
            17im//2835, -107//14175, -172im//155925, -109//133650]
        e22 = T[1, -im, -1, im//3, -1//3, 2im//15, 7//45, -2im//63, 8//315,
            -17im//2835, -107//14175, 172im//155925, -109//133650]
        e12 = T[0, -1, im, 1//3, -im//3, 2//15, -2im//15, -2//63, 2im//63,
            -17//2835, 17im//2835, 172//155925, -172im//155925]
        e21 = T[0, 1, im, -1//3, -im//3, -2//15, -2im//15, 2//63, 2im//63,
            17//2835, 17im//2835, -172//155925, -172im//155925]

        for k in 0:12
            @test coeffs[k+1] == [e11[k+1] e12[k+1]; e21[k+1] e22[k+1]]
        end

        # M(x)^2 = -x^2 (2 + x^2) I, which is what gives the closed form
        # exp(M) = cos(t) I + sin(t)/t M with t = x sqrt(2 + x^2).
        @test A^2 == 2 * one(A)
        @test B^2 == one(B)
        @test A * B == -B * A
    end

    @testset "order condition holds exactly for (M, N) = (k, k)" begin
        for k in 1:6
            Pc, Qc = matrix_pade_right_coeffs(coeffs, k, k)
            for j in 0:(2k)
                @test iszero(_torder_coeff(coeffs, Pc, Qc, j, :right))
            end
        end
    end

    @testset "exact rational evaluation at x = 1//2" begin
        expected = [
            (4 + 2im) // 5,
            (224 + 138im) // 305,
            (46864 + 29082im) // 64025,
            (268089248 + 166499290im) // 366396473,
            (21949205668979 + 13631856216360im) // 29997991820821,
            (6608375858586849397 + 4104226355542786536im) // 9031674702833386045,
        ]
        for k in 1:6
            g = matrix_pade_right(coeffs, k, k, 1 // 2)
            @test eltype(g) <: Complex{<:Rational}
            @test g[1, 1] == expected[k]
            # The series is symmetric under transpose-and-conjugate-i here.
            @test g[2, 2] == conj(expected[k])
        end
    end

    @testset "64-bit rationals overflow, BigInt rationals do not" begin
        small = _texp_iab_series(Complex{Rational{Int}}[1 im; -im -1],
            Complex{Rational{Int}}[0 1; 1 0], 12)
        # (1,1) and (2,2) still fit...
        @test matrix_pade_right(small, 2, 2, 1 // 2)[1, 1] == (224 + 138im) // 305
        # ...but (4,4) does not.
        @test_throws OverflowError matrix_pade_right(small, 4, 4, 1 // 2)
    end

    @testset "float accuracy improves with k and beats truncation" begin
        cf = [ComplexF64.(c) for c in coeffs]
        Af = ComplexF64.(A)
        Bf = ComplexF64.(B)
        exact(x) = exp(im * Af * x + im * Bf * x^2)

        for x in (0.5, 1.0)
            e = exact(x)
            errs = [norm(matrix_pade_right(cf, k, k, x) - e) for k in 1:6]
            @test issorted(errs, rev=true)
            # left and right agree closely on this series
            for k in 1:6
                @test norm(matrix_pade_left(cf, k, k, x) - e) ≈ errs[k] rtol = 1e-6
            end
            # from k = 4 on, the approximant is clearly ahead of truncating the
            # series at the same number of coefficients
            for k in 4:6
                trunc = sum(cf[j+1] * x^j for j in 0:(2k))
                @test errs[k] < norm(trunc - e)
            end
        end
    end

    @testset "observed order is x^(M+N+1)" begin
        cf = [ComplexF64.(c) for c in coeffs]
        Af = ComplexF64.(A)
        Bf = ComplexF64.(B)
        exact(x) = exp(im * Af * x + im * Bf * x^2)
        # Stop at k = 4: beyond that the float errors reach the roundoff floor
        # and the ratio stops being meaningful.
        for k in 1:4
            e1 = norm(matrix_pade_right(cf, k, k, 0.1) - exact(0.1))
            e2 = norm(matrix_pade_right(cf, k, k, 0.05) - exact(0.05))
            # The remainder is only asymptotically a pure power, so bracket the
            # ratio rather than pinning it (same convention as the other
            # order-decay tests here). Measured: 7.97, 32.04, 127.83, 510.37.
            @test 2.0^(2k) < e1 / e2 < 2.0^(2k + 2)
        end
    end
end
