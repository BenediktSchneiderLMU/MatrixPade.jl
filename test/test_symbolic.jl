using Test
using LinearAlgebra
using MatrixPade
using Symbolics

isdefined(@__MODULE__, :_texp_iab_series) || include("exp_series_testutils.jl")

# Both tracks are generic over the number type, so they work with Symbolics.jl
# either in the evaluation point or in the series coefficients themselves.
#
# These tests substitute a rational point back into the symbolic result and
# compare against direct exact evaluation, rather than asserting a particular
# `simplify` normal form, which is not stable across Symbolics versions.
@testset "Symbolics.jl interoperability" begin
    @variables z

    _sub(expr, val) = Symbolics.value(substitute(expr, Dict(z => val)))

    @testset "symbolic evaluation point, MPTA" begin
        # Gu (2004), Example 4.4: f(z) = exp(A z), A = [0 1; 0 -2].
        coeffs = [
            [1//1 0//1; 0//1 1//1],
            [0//1 1//1; 0//1 -2//1],
            [0//1 -1//1; 0//1 2//1],
            [0//1 2//3; 0//1 -4//3],
            [0//1 -1//3; 0//1 2//3],
            [0//1 2//15; 0//1 -4//15],
        ]
        sym = mpta(coeffs, 3, 2)(z)
        @test eltype(sym) <: Num

        # The (1,1) entry of exp(A z) is identically 1 and the (2,1) entry is
        # identically 0; the approximant reproduces both exactly, not just to
        # some order.
        @test isequal(Symbolics.value(simplify(sym[1, 1])), 1)
        @test iszero(Symbolics.value(simplify(sym[2, 1])))

        for pt in (1//2, 3//1, -2//7)
            @test _sub.(sym, pt) == mpta(coeffs, 3, 2, pt)
        end
    end

    @testset "symbolic evaluation point, matrix Pade form" begin
        # Beckermann and Labahn (1994), section 5.
        coeffs = [
            [1 0; 0 1], [0 0; 0 0], [1 0; 0 1],
            [0 0; 0 0], [2 0; 0 1], [-1 0; -1 0],
        ]
        sym = matrix_pade_left(coeffs, 2, 3)(z)
        @test eltype(sym) <: Num

        for pt in (1//2, 3//1, -2//7)
            @test _sub.(sym, pt) == matrix_pade_left(coeffs, 2, 3, pt)
        end
    end

    @testset "symbolic evaluation point, complex series" begin
        T = Complex{Rational{BigInt}}
        coeffs = _texp_iab_series(_treadme_A(T), _treadme_B(T), 12)
        sym = matrix_pade_right(coeffs, 2, 2)(z)

        # Entries are Complex{Num}: substitute into the parts.
        subc(w, val) = _sub(real(w), val) + im * _sub(imag(w), val)
        @test subc.(sym, 1//2) == matrix_pade_right(coeffs, 2, 2, 1//2)
        @test subc(sym[1, 1], 1//2) == (224 + 138im)//305
    end

    @testset "symbolic series coefficients" begin
        @variables a
        coeffs = [[1 0; 0 1], [0 a; 0 0], [0 0; a 0]]
        qcoeffs, Pcoeffs = mpta_coeffs(coeffs, 1, 1)
        @test isequal(Symbolics.value(simplify(qcoeffs[1])), Symbolics.value(a^2))
        @test iszero(Symbolics.value(simplify(qcoeffs[2])))
        @test length(Pcoeffs) == 2
    end
end
