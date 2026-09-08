module MatrixPade

using LinearAlgebra

export matdot, mpta, mpta_coeffs, PadeApproximant
export fphps, defect, matrix_pade, matrix_pade_coeffs, matrix_pade_right,
    matrix_pade_right_coeffs, matrix_pade_left, matrix_pade_left_coeffs,
    MatrixPadeForm

include("scalar_product.jl")
include("mpta.jl")
include("fphps.jl")
include("matrix_pade.jl")

end # module MatrixPade
