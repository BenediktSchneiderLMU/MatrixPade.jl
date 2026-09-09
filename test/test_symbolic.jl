using Test
using LinearAlgebra
using MatrixPade
using Symbolics

isdefined(@__MODULE__, :_texp_iab_series) || include("exp_series_testutils.jl")

# The matrix Pade forms are generic over the number type, so they accept a
# Symbolics.jl evaluation point and hand back the rational function itself.
#
# These tests substitute a rational point back into the symbolic result and
# compare against direct exact evaluation, rather than asserting a particular
# `simplify` normal form, which is not stable across Symbolics versions.
@testset "Symbolics.jl interoperability" begin
    @variables z

    _sub(expr, val) = Symbolics.value(substitute(expr, Dict(z => val)))

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
end
