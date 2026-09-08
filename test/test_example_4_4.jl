using Test
using LinearAlgebra
using MatrixPade

# Example 4.4 of Gu, "Matrix Padé-type approximant and directional matrix
# Padé approximant in the inner product space", J. Comput. Appl. Math.
# 164-165 (2004) 365-385. Seeks (3/2)_f(z) for f(z) = e^{Az},
# A = [0 1; 0 -2].
@testset "Example 4.4 (Gu 2004)" begin
    T = Rational{Int}
    c0 = T[1 0; 0 1]
    c1 = T[0 1; 0 -2]
    c2 = T[0 -1; 0 2]
    c3 = T[0 2//3; 0 -4//3]
    c4 = T[0 -1//3; 0 2//3]
    c5 = T[0 2//15; 0 -4//15]
    coeffs = [c0, c1, c2, c3, c4, c5]

    m, n = 3, 2

    # Sanity check, computed independently of the package's internals:
    # reproduce the paper's stated det{H2(c2)} = 50/27.
    H = [matdot(c2, c2) matdot(c2, c3)
         matdot(c3, c3) matdot(c3, c4)]
    @test det(H) == 50 // 27

    qcoeffs, Pcoeffs = mpta_coeffs(coeffs, m, n)

    @test qcoeffs == T[50//27, 40//27, 10//27]

    expected_P = [
        T[50//27 0; 0 50//27],
        T[40//27 50//27; 0 -20//9],
        T[10//27 -10//27; 0 10//9],
        T[0 10//81; 0 -20//81],
    ]
    @test Pcoeffs == expected_P

    # Independent residual check: P must equal the degree-m truncation of
    # q(z)*f(z), i.e. P_j = sum_{k=0}^{min(j,n)} q_k * c_{j-k}. This does
    # not rely on the hand-derived expected values above.
    for j in 0:m
        conv = zero(c0)
        for k in 0:min(j, n)
            conv = conv + qcoeffs[k+1] * coeffs[j-k+1]
        end
        @test conv == Pcoeffs[j+1]
    end
end
