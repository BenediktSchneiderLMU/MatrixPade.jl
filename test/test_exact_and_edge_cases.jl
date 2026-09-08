using Test
using MatrixPade

@testset "exactness and edge cases" begin
    # A small, non-degenerate (0/1) example: c0 = I, c1 = 2I.
    c0i = [1 0; 0 1]
    c1i = [2 0; 0 2]
    coeffs_int = [c0i, c1i]

    @testset "integer input promotes to exact Rational" begin
        qcoeffs, Pcoeffs = mpta_coeffs(coeffs_int, 0, 1)
        @test eltype(qcoeffs) <: Rational
        @test eltype(eltype(Pcoeffs)) <: Rational
        @test qcoeffs == [2 // 1, -4 // 1]
        @test Pcoeffs == [[2 // 1 0; 0 2 // 1]]
    end

    @testset "float input stays float" begin
        c0f = Float64[1 0; 0 1]
        c1f = Float64[2 0; 0 2]
        coeffs_float = [c0f, c1f]
        qcoeffs_f, Pcoeffs_f = mpta_coeffs(coeffs_float, 0, 1)
        @test eltype(qcoeffs_f) <: AbstractFloat
        @test eltype(eltype(Pcoeffs_f)) <: AbstractFloat
    end

    @testset "validation" begin
        # needs m+n+1 = 6 coefficients, only 2 given
        @test_throws ArgumentError mpta_coeffs(coeffs_int, 3, 2)
        # n must be >= 1
        @test_throws ArgumentError mpta_coeffs(coeffs_int, 1, 0)
        # m must be >= n-1
        @test_throws ArgumentError mpta_coeffs(coeffs_int, 0, 2)
    end
end
