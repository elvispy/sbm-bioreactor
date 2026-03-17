using Test
using SBM_Bioreactor

@testset "Module Loading" begin
    @test true
end

include("test_navier_stokes.jl")
