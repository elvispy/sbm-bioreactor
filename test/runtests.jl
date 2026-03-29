# Main test runner for the SBM_Bioreactor package.
# Fast checks run by default; slow scientific checks are opt-in.
using Test
using SBM_Bioreactor

@testset "Module Loading" begin
    # Verify the package exports and module environment load correctly.
    @test true
end

# Fast default suite: keep this aligned with the current explicit-BDF1 and blocked-solver path.
# Each included file contains dedicated @testset blocks for specific components.
include("test_navier_stokes.jl")
include("test_rheology.jl")
include("test_migration.jl")
include("test_solver_fast.jl")
include("test_examples.jl")

# Slow scientific checks are opt-in so the default suite stays responsive.
if get(ENV, "SBM_RUN_SLOW_TESTS", "0") == "1"
    include("test_solver.jl")
    include("test_mms.jl")
    include("test_mms_unsteady.jl")
end
