using Test
using LinearAlgebra
using MatrixPade

# Shared, independent polynomial helpers; guarded so that running several of
# these files in one session does not redefine them.
isdefined(@__MODULE__, :_torder_coeff) || include("matrix_pade_testutils.jl")

# Section 5 of Beckermann and Labahn (1994), continued: the (2, 3) left-hand
# matrix Pade forms of the same A(z), p. 816-817.
@testset "Section 5, left-hand matrix Pade form" begin
    T = Rational{Int}
    coeffs = [
        T[1 0; 0 1],    # A_0
        T[0 0; 0 0],    # A_1
        T[1 0; 0 1],    # A_2
        T[0 0; 0 0],    # A_3
        T[2 0; 0 1],    # A_4
        T[-1 0; -1 0],  # A_5
    ]
    M, N = 2, 3

    # Table 1, row 2.2: s = q = 2, n = (M, M, N, N), sigma = q(M+N+1) = 12,
    # and F as printed on p. 816. Built independently of the package.
    sigma = 2 * (M + N + 1)
    f1 = zeros(T, sigma); f1[1] = 1                        # 1
    f2 = zeros(T, sigma); f2[2] = 1                        # z
    f3 = zeros(T, sigma)                                   # -1-z^4-2z^8+z^10
    f3[1] = -1; f3[5] = -1; f3[9] = -2; f3[11] = 1
    f4 = zeros(T, sigma)                                   # -z-z^5-z^9+z^10
    f4[2] = -1; f4[6] = -1; f4[10] = -1; f4[11] = 1

    F = [f1, f2, f3, f4]
    n = [M, M, N, N]
    B, d, _ = fphps(F, sigma, n, 2)

    # "In this case the defects are 1, -1, 1, and 1, respectively."
    @test defect(d) == [1, -1, 1, 1]
    for l in 1:4
        row = [B[l, lam] for lam in 1:4]
        @test all(iszero, _tresidual(row, F, 2, sigma))  # ord P_l >= sigma
        @test _tdefect(row, n) == defect(d)[l]           # dct P_l = d_l + 1
    end
    # "so the solution space L_0^12 is of the form a P_1 + b P_3 + c P_4":
    # dimension 3, while only r = p = 2 rows are needed.
    @test sum(max(dl + 1, 0) for dl in d) == 3

    Pcoeffs, Qcoeffs = matrix_pade_left_coeffs(coeffs, M, N)

    @test length(Pcoeffs) == M + 1
    @test length(Qcoeffs) == N + 1

    # Our selection rule takes the a- and b-generators in canonical order, so
    #   P(z) = [-1+z^2  1 ;  -z      -1     ]
    #   Q(z) = [-1+2z^2 1-z^2; -z+z^3 -1+z^2]
    # which is the paper's displayed form with its two rows interchanged (the
    # form is not unique here -- see the paper's own remark on p. 817).
    @test Pcoeffs == [T[-1 1; 0 -1], T[0 0; -1 0], T[1 0; 0 0]]
    @test Qcoeffs == [T[-1 1; 0 -1], T[0 0; -1 0], T[2 -1; 0 1], T[0 0; 1 0]]

    # Independent check of the defining relation Q(z) A(z) - P(z) = O(z^6).
    for j in 0:(M+N)
        @test iszero(_torder_coeff(coeffs, Pcoeffs, Qcoeffs, j, :left))
    end
    @test any(!iszero, Qcoeffs)

    g = matrix_pade_left(coeffs, M, N)
    @test g.side === :left
    @test g.P == Pcoeffs
    @test g.Q == Qcoeffs

    # "Note that the denominator has a nonzero determinant, indeed that Q(0)
    # is nonsingular. Therefore, unlike the case for approximants on the
    # right, one can always form the rational expression Q(z)^-1 P(z)."
    @test det(Qcoeffs[1]) != 0
    @test g(0) == coeffs[1]

    # Evaluation stays exact, and by construction Q(z) g(z) == P(z).
    for z in (1 // 2, 3, -2 // 7)
        gz = g(z)
        @test eltype(gz) <: Rational
        Qz = sum(Qcoeffs[k+1] * z^k for k in 0:N)
        Pz = sum(Pcoeffs[k+1] * z^k for k in 0:M)
        @test Qz * gz == Pz
    end

    # The evaluated error really is O(z^{M+N+1}): halving z shrinks it by a
    # factor tending to 2^(M+N+1) = 64, approached from below, so test a band.
    Atrunc(z) = sum(coeffs[k+1] * z^k for k in 0:(M + N))
    err(z) = maximum(abs, g(z) - Atrunc(z))
    @test err(1 // 10) > 0
    @test 2^(M + N) < err(1 // 10) / err(1 // 20) < 2^(M + N + 2)
end
