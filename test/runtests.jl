using Test
using LinearAlgebra
using MatrixPade

@testset "MatrixPade.jl" begin
    @testset "matrix_pade_eval" begin
        coeffs = [Matrix{Float64}(I, 2, 2), zeros(2, 2)]
        @test matrix_pade_eval(coeffs, 0.5) == Matrix{Float64}(I, 2, 2)
        @test_throws ArgumentError matrix_pade_eval(Matrix{Float64}[], 0.5)
    end
end
