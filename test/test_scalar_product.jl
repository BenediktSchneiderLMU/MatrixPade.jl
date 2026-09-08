using Test
using LinearAlgebra
using MatrixPade

@testset "matdot" begin
    A = [1 2; 3 4]
    B = [5 6; 7 8]

    @test matdot(A, B) == sum(A .* B)
    @test matdot(A, B) == 1 * 5 + 2 * 6 + 3 * 7 + 4 * 8

    # bilinearity
    A2 = [0 1; 1 0]
    @test matdot(A + A2, B) == matdot(A, B) + matdot(A2, B)
    k = 3
    @test matdot(k .* A, B) == k * matdot(A, B)

    # symmetry: (A,B) == (B,A)
    @test matdot(A, B) == matdot(B, A)

    # no conjugation for complex entries (this is the paper's departure
    # from the usual Hermitian inner product)
    C = ComplexF64[1+2im 0; 0 0]
    @test matdot(C, C) == sum(C .^ 2)
    @test matdot(C, C) != norm(C)^2

    # exact on Rationals
    Aq = Rational{Int}[1//2 1//3; 1//4 1//5]
    Bq = Rational{Int}[2//1 3//1; 4//1 5//1]
    @test matdot(Aq, Bq) == sum(Aq .* Bq)
    @test matdot(Aq, Bq) isa Rational
end
