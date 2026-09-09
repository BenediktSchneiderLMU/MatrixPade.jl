using Test
using LinearAlgebra
using MatrixPade

@testset "MatrixPade.jl" begin
    include("test_fphps.jl")
    include("test_matrix_pade_right.jl")
    include("test_matrix_pade_left.jl")
    include("test_matrix_pade_types.jl")
    include("test_readme_example.jl")
    include("test_symbolic.jl")
end
