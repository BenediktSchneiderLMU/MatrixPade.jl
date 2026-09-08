using Test
using LinearAlgebra
using MatrixPade

@testset "MatrixPade.jl" begin
    include("test_scalar_product.jl")
    include("test_example_4_4.jl")
    include("test_pade_approximant.jl")
    include("test_exact_and_edge_cases.jl")
    include("test_fphps.jl")
    include("test_matrix_pade_right.jl")
    include("test_matrix_pade_left.jl")
    include("test_matrix_pade_types.jl")
    include("test_mpta_degenerate.jl")
    include("test_readme_example.jl")
    include("test_symbolic.jl")
end
