using Test
using MatrixPade

# Shared, independent polynomial helpers; guarded so that running several of
# these files in one session does not redefine them.
isdefined(@__MODULE__, :_torder_coeff) || include("matrix_pade_testutils.jl")

# Example 4.4 of Beckermann and Labahn, "A uniform approach for the fast
# computation of matrix-type Pade approximants", SIAM J. Matrix Anal. Appl.
# 15 (1994) 804-823, p. 813. Ground truth for the scalar solver itself,
# independent of the matrix Pade layer:
#
#   m = 4, s = 2, n = (2,2,2,2),
#   F = (1, z, z/(1-z^4) + z^10, z/(1+z^4) + z^12)^T + O(z^16).
@testset "FPHPS (Beckermann-Labahn Example 4.4)" begin
    T = Rational{BigInt}
    L = 16

    f1 = zeros(T, L); f1[1] = 1
    f2 = zeros(T, L); f2[2] = 1
    f3 = zeros(T, L)
    f4 = zeros(T, L)
    for k in 0:L
        idx = 4k + 2          # z^(4k+1)
        idx <= L || break
        f3[idx] += 1          # z/(1-z^4)
        f4[idx] += (-1)^k     # z/(1+z^4)
    end
    f3[11] += 1               # + z^10
    f4[13] += 1               # + z^12

    F = [f1, f2, f3, f4]
    n = [2, 2, 2, 2]

    # "An application of FPHPS gives the values pi_0, pi_1, ..., pi_13 =
    #  1, 2, 1, 2, 1, 3, 1, 3, 1, 4, 2, 4, 3, 4."
    _, _, pivots = fphps(F, 14, n, 2)
    @test pivots == [1, 2, 1, 2, 1, 3, 1, 3, 1, 4, 2, 4, 3, 4]

    P, d, _ = fphps(F, 10, n, 2)

    # The sigma-basis for sigma = 10 printed on p. 813, rows P_{l,10}.
    @test [_ttrim(P[1, lam]) for lam in 1:4] == [T[0, 0, 0, 0, 0, 1], T[], T[], T[]]
    @test [_ttrim(P[2, lam]) for lam in 1:4] == [T[], T[0, 0, 1], T[-1//2], T[1//2]]
    @test [_ttrim(P[3, lam]) for lam in 1:4] == [T[], T[1, 0, -1], T[-1//2, 0, 1], T[-1//2]]
    @test [_ttrim(P[4, lam]) for lam in 1:4] == [T[], T[0, -2], T[0, 1], T[0, 1]]

    # "The defects for this basis are -2, 1, 1 and 2, respectively."
    @test defect(d) == [-2, 1, 1, 2]
    for l in 1:4
        @test _tdefect([P[l, lam] for lam in 1:4], n) == defect(d)[l]
    end

    # The s-residuals P_10(z^2) . F(z) printed alongside them, checked as far
    # as the paper states them explicitly.
    R = [_tresidual([P[l, lam] for lam in 1:4], F, 2, L) for l in 1:4]
    @test all(iszero, R[l][1:10] for l in 1:4)          # every row has order >= 10
    @test R[1][11:16] == T[1, 0, 0, 0, 0, 0]            # z^10 + O(z^26)
    @test R[2][11:14] == T[-1//2, 0, 1//2, -1]          # -z^10/2 + z^12/2 - z^13
    @test R[3][11:15] == T[-1//2, 0, -1//2, 1, 1]       # -z^10/2 - z^12/2 + z^13 + z^14
    @test R[4][11:16] == T[0, 2, 1, 0, 1, 0]            # 2z^11 + z^12 + z^14

    # sigma = 0 returns the initial basis unchanged.
    P0, d0, piv0 = fphps(F, 0, n, 2)
    @test d0 == n
    @test isempty(piv0)
    @test all(_ttrim(P0[l, lam]) == (l == lam ? T[1] : T[]) for l in 1:4, lam in 1:4)

    # Integer input is promoted so the divisions stay exact.
    Pi, _, _ = fphps([Int.(f) for f in F], 10, n, 2)
    @test eltype(Pi[2, 3]) <: Rational

    @testset "validation" begin
        @test_throws ArgumentError fphps([f1], 4, [2], 2)          # m < 2
        @test_throws ArgumentError fphps(F, 4, n, 0)               # s < 1
        @test_throws ArgumentError fphps(F, -1, n, 2)              # sigma < 0
        @test_throws ArgumentError fphps(F, 4, [2, 2, 2], 2)       # length(n) != m
        @test_throws ArgumentError fphps(F, 20, n, 2)              # F too short
    end
end
