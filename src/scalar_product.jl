"""
    matdot(A, B)

Bilinear scalar product of two same-shape matrices, `sum(A .* B)`. Unlike
the usual Hermitian inner product for complex matrices, this does **not**
conjugate either argument — this is the scalar product Gu (2004) defines
for the construction of the matrix Padé-type approximant.
"""
matdot(A::AbstractMatrix, B::AbstractMatrix) = sum(A .* B)
