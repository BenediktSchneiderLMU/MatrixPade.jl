using Test
using LinearAlgebra
using MatrixPade

# Shared, independent polynomial helpers; guarded so that running several of
# these files in one session does not redefine them.
isdefined(@__MODULE__, :_torder_coeff) || include("matrix_pade_testutils.jl")

@testset "matrix Pade forms: types, shapes and edge cases" begin
    # A wide (2 by 3) series suits the right-hand form, a tall (3 by 2) one
    # the left-hand form: dim L_0^sigma >= q(N+1) - pN on the right and
    # p(N+1) - qN on the left.
    wide = [[1 0 2; 0 1 1], [0 1 0; 1 0 -1], [1 1 0; 0 2 1]]
    tall = [[1 0; 0 1; 2 1], [0 1; 1 0; 0 -1], [1 1; 0 2; 1 0]]
    M, N = 1, 1

    @testset "integer input promotes to exact Rational" begin
        for (A, side) in ((wide, :right), (tall, :left))
            Pc, Qc = matrix_pade_coeffs(A, M, N; side=side)
            @test eltype(Pc[1]) <: Rational
            @test eltype(Qc[1]) <: Rational
            @test eltype(matrix_pade(A, M, N, 1 // 3; side=side)) <: Rational
            # Integer evaluation points stay exact too.
            @test eltype(matrix_pade(A, M, N, 2; side=side)) <: Rational
        end
    end

    @testset "float input stays float" begin
        Af = [Float64.(c) for c in wide]
        Pc, Qc = matrix_pade_right_coeffs(Af, M, N)
        @test eltype(Pc[1]) === Float64
        @test eltype(Qc[1]) === Float64
        @test eltype(matrix_pade_right(Af, M, N, 0.25)) === Float64
    end

    @testset "square denominator (default r)" begin
        # right: P is p by q, Q is q by q; left: P is p by q, Q is p by p.
        for (A, side, pshape, qshape) in
            ((wide, :right, (2, 3), (3, 3)), (tall, :left, (3, 2), (3, 3)))
            Pc, Qc = matrix_pade_coeffs(A, M, N; side=side)
            @test size(Pc[1]) == pshape
            @test size(Qc[1]) == qshape
            @test length(Pc) == M + 1 && length(Qc) == N + 1
            @test any(!iszero, Qc)
            for j in 0:(M+N)
                @test iszero(_torder_coeff(A, Pc, Qc, j, side))
            end
            # A square Q makes the fraction reproduce A to order M+N+1.
            f = matrix_pade(A, M, N; side=side)
            @test size(f(1 // 3)) == size(A[1])
            Atrunc(z) = sum(A[k+1] * z^k for k in 0:(M+N))
            err(z) = maximum(abs, f(z) - Atrunc(z))
            @test err(1 // 10) > 0
            @test 2^(M + N) < err(1 // 10) / err(1 // 20) < 2^(M + N + 2)
        end
    end

    @testset "rectangular denominator uses the pseudoinverse" begin
        # r below the default leaves Q rectangular: tall on the right,
        # wide on the left. The order condition still holds; the fraction
        # no longer reproduces A, but must still evaluate to A's shape.
        for (A, side, r, qshape) in
            ((wide, :right, 2, (3, 2)), (tall, :left, 2, (2, 3)))
            Pc, Qc = matrix_pade_coeffs(A, M, N; side=side, r=r)
            @test size(Qc[1]) == qshape
            @test any(!iszero, Qc)
            for j in 0:(M+N)
                @test iszero(_torder_coeff(A, Pc, Qc, j, side))
            end
            f = matrix_pade(A, M, N; side=side, r=r)
            v = f(1 // 3)
            @test size(v) == size(A[1])
            @test eltype(v) <: Rational   # pinv stays exact, unlike LinearAlgebra.pinv
        end
    end

    @testset "wrappers agree with the keyword form" begin
        @test matrix_pade_right_coeffs(wide, M, N) ==
              matrix_pade_coeffs(wide, M, N; side=:right)
        @test matrix_pade_left_coeffs(tall, M, N) ==
              matrix_pade_coeffs(tall, M, N; side=:left)
        f = matrix_pade_right(wide, M, N)
        @test f.P == matrix_pade_right_coeffs(wide, M, N)[1]
        @test f.Q == matrix_pade_right_coeffs(wide, M, N)[2]
        @test matrix_pade_right(wide, M, N, 1 // 3) == f(1 // 3)
        @test matrix_pade_left(tall, M, N, 1 // 3) == matrix_pade_left(tall, M, N)(1 // 3)
    end

    @testset "validation" begin
        @test_throws ArgumentError matrix_pade_coeffs(wide, M, N; side=:up)   # bad side
        @test_throws ArgumentError matrix_pade_coeffs(wide, -1, N)            # M < 0
        @test_throws ArgumentError matrix_pade_coeffs(wide, M, -1)            # N < 0
        @test_throws ArgumentError matrix_pade_coeffs(wide, 3, 3)             # too few coeffs
        @test_throws ArgumentError matrix_pade_coeffs(wide, M, N; r=0)        # r < 1
        # A wide series has no square left-hand denominator here (p < q).
        @test_throws ArgumentError matrix_pade_coeffs(wide, M, N; side=:left)
        # ...and the message says so, rather than failing obscurely.
        err = try
            matrix_pade_coeffs(wide, M, N; side=:left)
        catch e
            e
        end
        @test occursin("p >= q", err.msg)
        # Mismatched coefficient sizes.
        @test_throws ArgumentError matrix_pade_coeffs([[1 0; 0 1], [1 0 0; 0 1 0]], 0, 0)
    end
end
