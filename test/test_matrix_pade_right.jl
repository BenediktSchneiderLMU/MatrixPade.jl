using Test
using MatrixPade

# Shared, independent polynomial helpers; guarded so that running several of
# these files in one session does not redefine them.
isdefined(@__MODULE__, :_torder_coeff) || include("matrix_pade_testutils.jl")

# Section 5 of Beckermann and Labahn, "A uniform approach for the fast
# computation of matrix-type Pade approximants", SIAM J. Matrix Anal. Appl.
# 15 (1994) 804-823. Seeks the (2, 3) right-hand matrix Pade form of
#
#   A(z) = [1+z^2+2z^4-z^5+z^6   z^7                ] + O(z^8).
#          [-z^5                 1+z^2+z^4+z^7      ]
@testset "Section 5, right-hand matrix Pade form" begin
    T = Rational{Int}
    # Only A_0 ... A_{M+N} = A_0 ... A_5 are needed.
    coeffs = [
        T[1 0; 0 1],    # A_0
        T[0 0; 0 0],    # A_1
        T[1 0; 0 1],    # A_2
        T[0 0; 0 0],    # A_3
        T[2 0; 0 1],    # A_4
        T[-1 0; -1 0],  # A_5
    ]
    M, N = 2, 3

    # The scalar problem the paper reduces this to (Table 1, row 2.1):
    # s = p = 2, n = (M, M, N, N), sigma = p(M+N+1) = 12, and F as printed
    # on p. 815. Built here independently of the package's internals.
    sigma = 2 * (M + N + 1)
    f1 = zeros(T, sigma); f1[1] = 1                       # 1
    f2 = zeros(T, sigma); f2[2] = 1                       # z
    # -1 - z^4 - 2z^8 + z^10 + z^11 - z^12
    f3 = zeros(T, sigma)
    f3[1] = -1; f3[5] = -1; f3[9] = -2; f3[11] = 1; f3[12] = 1
    # -z - z^5 - z^9 - z^14 - z^15
    f4 = zeros(T, sigma)
    f4[2] = -1; f4[6] = -1; f4[10] = -1

    F = [f1, f2, f3, f4]
    n = [M, M, N, N]
    B, d, _ = fphps(F, sigma, n, 2)

    # A sigma-basis is not unique: Definition 3.2 pins the solution space,
    # not a representative. (Ours differs from the one printed on p. 815 by a
    # single elementary row operation -- our row 3 minus our row 1 is exactly
    # their row 3.) So check the properties the paper states, not the printed
    # matrix.
    #
    # "The defects for this basis are 0, 0, 0, and 2, respectively."
    @test defect(d) == [0, 0, 0, 2]
    for l in 1:4
        row = [B[l, lam] for lam in 1:4]
        @test all(iszero, _tresidual(row, F, 2, sigma))  # ord P_l >= sigma
        @test _tdefect(row, n) == defect(d)[l]           # dct P_l = d_l + 1
    end
    # dim L_0^12 = sum_l max(d_l + 1, 0) = 2, and the paper's generator of
    # that space is P_{4,12} = [0, -1, 0, -1 + z^2].
    @test sum(max(dl + 1, 0) for dl in d) == 2
    @test [B[4, lam] for lam in 1:4] == [T[0, 0, 0], T[-1, 0, 0], T[0, 0, 0], T[-1, 0, 1]]

    Pcoeffs, Qcoeffs = matrix_pade_right_coeffs(coeffs, M, N)

    @test length(Pcoeffs) == M + 1
    @test length(Qcoeffs) == N + 1

    # P(z) = [0 0; -1 -z],  Q(z) = [0 0; -1+z^2  -z+z^3]
    @test Pcoeffs == [T[0 0; -1 0], T[0 0; 0 -1], T[0 0; 0 0]]
    @test Qcoeffs == [T[0 0; -1 0], T[0 0; 0 -1], T[0 0; 1 0], T[0 0; 0 1]]

    # Independent check of the defining relation, computed here by plain
    # truncated polynomial arithmetic rather than through the package:
    # A(z) Q(z) - P(z) must vanish to order M+N+1 = 6.
    for j in 0:(M+N)
        @test iszero(_torder_coeff(coeffs, Pcoeffs, Qcoeffs, j, :right))
    end
    # ...and it is not vanishing merely because P and Q are trivial.
    @test any(!iszero, Qcoeffs)

    # The paper notes that this Q(z) is singular for every z, so no right
    # matrix Pade fraction of type (2, 3) exists here.
    f = matrix_pade_right(coeffs, M, N)
    @test f.side === :right
    @test_throws ArgumentError f(1 // 2)

    # Byproducts listed on p. 817: the right matrix Pade forms of type
    # (0, 1) and (1, 2) satisfy their own order conditions.
    for (Mi, Ni) in ((0, 1), (1, 2))
        Pi, Qi = matrix_pade_right_coeffs(coeffs, Mi, Ni)
        @test length(Pi) == Mi + 1
        @test length(Qi) == Ni + 1
        @test any(!iszero, Qi)
        for j in 0:(Mi+Ni)
            @test iszero(_torder_coeff(coeffs, Pi, Qi, j, :right))
        end
    end
end
