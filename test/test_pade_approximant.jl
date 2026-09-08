using Test
using MatrixPade

@testset "PadeApproximant / mpta" begin
    T = Rational{Int}
    c0 = T[1 0; 0 1]
    c1 = T[0 1; 0 -2]
    c2 = T[0 -1; 0 2]
    c3 = T[0 2//3; 0 -4//3]
    c4 = T[0 -1//3; 0 2//3]
    c5 = T[0 2//15; 0 -4//15]
    coeffs = [c0, c1, c2, c3, c4, c5]
    m, n = 3, 2

    qcoeffs, Pcoeffs = mpta_coeffs(coeffs, m, n)
    pa = mpta(coeffs, m, n)

    @test pa isa PadeApproximant
    @test pa.q == qcoeffs
    @test pa.P == Pcoeffs

    # Horner evaluation, written independently of the package internals.
    qz(coefs, z) = sum(coefs[i+1] * z^i for i in 0:length(coefs)-1)
    Pz(coefs, z) = sum(coefs[i+1] * z^i for i in 0:length(coefs)-1)

    z = 1 // 3
    expected = Pz(Pcoeffs, z) ./ qz(qcoeffs, z)

    @test pa(z) == expected
    @test mpta(coeffs, m, n, z) == expected
end
