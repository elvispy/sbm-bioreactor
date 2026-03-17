using Test
using SBM_Bioreactor

@testset "Module Loading" begin
    @test true
end

include("test_navier_stokes.jl")
include("test_rheology.jl")
include("test_migration.jl")
include("test_solver.jl")
