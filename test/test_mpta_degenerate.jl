using Test
using LinearAlgebra
using MatrixPade

isdefined(@__MODULE__, :_texp_iab_series) || include("exp_series_testutils.jl")

# Gu's scalar product is bilinear and does not conjugate, so it can be
# isotropic on a complex series. For f(x) = exp(i A x + i B x^2) with the
# README's A and B, every odd coefficient is a null vector of matdot, the
# Hankel matrix is singular, and no Pade-type denominator exists.
@testset "MPTA degeneracy on an isotropic series" begin
    T = Complex{Rational{BigInt}}
    coeffs = _texp_iab_series(_treadme_A(T), _treadme_B(T), 12)

    @testset "matdot vanishes whenever either index is odd" begin
        for i in 0:12, j in 0:12
            if isodd(i) || isodd(j)
                @test iszero(matdot(coeffs[i+1], coeffs[j+1]))
            end
        end
        # ...and does not vanish identically on the even ones.
        @test !iszero(matdot(coeffs[1], coeffs[1]))
        @test !iszero(matdot(coeffs[1], coeffs[3]))
    end

    @testset "diagonal (k/k) approximants all degenerate" begin
        for k in 1:6
            @test_throws ArgumentError mpta_coeffs(coeffs, k, k)
            @test_throws ArgumentError mpta(coeffs, k, k)
        end
    end

    # A scan over every valid (m, n) with m+n+1 <= 13 leaves only these four.
    @testset "the surviving (m, n) still work" begin
        for (m, n) in ((0, 1), (6, 1), (8, 1), (10, 1))
            qcoeffs, Pcoeffs = mpta_coeffs(coeffs, m, n)
            @test !all(iszero, qcoeffs)
            @test length(qcoeffs) == n + 1
            @test length(Pcoeffs) == m + 1
        end
    end

    @testset "the matrix-denominator forms are unaffected" begin
        # Same series, genuine matrix denominator: no degeneracy.
        for k in 1:6
            f = matrix_pade_right(coeffs, k, k)
            @test f(0) == coeffs[1]
        end
    end
end
