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

    # q(z) f(z) - P(z) = O(z^{m+1}) is a guarantee about the numerator degree
    # alone. On a series with no special structure the bound is attained
    # exactly: the first nonvanishing residual coefficient sits at m+1 whatever
    # n is. (Gu's own exp(Az) example reaches m+n+1 instead, but only because
    # A = [0 1; 0 -2] satisfies A^2 = -2A and is rank-1 degenerate.)
    @testset "the O(z^{m+1}) bound is sharp" begin
        generic = [
            [-3 3; 2 0], [5 4; -5 4], [-1 4; -3 2], [4 3; 2 5], [3 3; 5 -2],
            [-2 1; 5 5], [-4 -4; -4 2], [-2 1; -5 1], [1 -3; 2 -2], [-3 5; -4 -3],
            [-5 5; -1 1], [-4 -3; 5 -2], [-5 0; 1 5], [4 1; 2 1], [5 0; -2 1],
            [-5 -4; -3 -5], [3 -4; 4 -3],
        ]
        c = [Rational{BigInt}.(x) for x in generic]

        # Index of the first nonvanishing coefficient of q(z) f(z) - P(z),
        # by plain truncated polynomial arithmetic.
        function residual_order(m, n)
            qcoeffs, Pcoeffs = mpta_coeffs(c, m, n)
            for j in 0:(length(c)-1)
                acc = zero(c[1])
                for k in 0:min(j, n)
                    acc = acc + qcoeffs[k+1] * c[j-k+1]
                end
                j <= m && (acc = acc - Pcoeffs[j+1])
                iszero(acc) || return j
            end
            return -1
        end

        for (m, n) in ((1, 1), (2, 1), (2, 2), (3, 2), (4, 2), (4, 3), (5, 3), (6, 4), (3, 4))
            @test residual_order(m, n) == m + 1
        end
    end
end
