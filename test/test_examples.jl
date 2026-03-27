using Test
using Gridap
using SBM_Bioreactor

@testset "HARV Example API" begin
    @test isdefined(SBM_Bioreactor, :build_harv_2d_case)
    @test isdefined(SBM_Bioreactor, :run_bioreactor_simulation)

    case = build_harv_2d_case(partition=(2, 2), dt=0.25, total_time=0.5)

    @test num_fields(case.X) == 5
    @test num_fields(case.Y) == 5
    @test case.metadata.nsteps == 2
    @test case.metadata.partition == (2, 2)
    @test case.metadata.dt == 0.25
    @test case.params.Φ0(Point(0.0, 0.03)) ≈ 0.1
    @test case.params.Φ0(Point(0.0, 0.0)) ≈ 0.0
    @test case.params.Γ0(Point(0.0, 0.0)) > 0.0

    blocked_case = build_harv_2d_case(partition=(2, 2), dt=0.25, total_time=0.5, blocked=true)
    blocked_state = interpolate_everywhere(
        [blocked_case.params.u0, blocked_case.params.p0, blocked_case.params.Φ0, blocked_case.params.C0, blocked_case.params.Γ0],
        blocked_case.X,
    )
    @test blocked_case.metadata.blocked == true
    @test num_fields(blocked_case.X) == 5
    @test occursin("BlockVector", string(typeof(get_free_dof_values(blocked_state))))

    result = run_bioreactor_simulation(
        case.X,
        case.Y,
        case.dΩ,
        case.metadata.dt,
        case.params,
        0;
        write_vtk_interval=0,
        collect_history=true,
    )

    @test length(result.history) == 1
    @test result.times == [0.0]
    @test length(result.final_state) == 5
    @test !haskey(result, :profile)

    result_no_history = run_bioreactor_simulation(
        case.X,
        case.Y,
        case.dΩ,
        case.metadata.dt,
        case.params,
        0;
        write_vtk_interval=0,
        collect_history=false,
    )

    @test num_fields(result_no_history) == 5

    profiled = run_bioreactor_simulation(
        case.X,
        case.Y,
        case.dΩ,
        case.metadata.dt,
        case.params,
        0;
        write_vtk_interval=0,
        collect_history=true,
        profile_steps=true,
    )

    @test haskey(profiled, :profile)
    @test profiled.profile.initial_setup_time >= 0.0
    @test isempty(profiled.profile.steps)
end

@testset "Visualization Helper Loading" begin
    plots_available = false
    try
        @eval using Plots
        plots_available = true
    catch
        plots_available = false
    end

    if plots_available
        include("../scripts/visualize.jl")
        case = build_harv_2d_case(partition=(2, 2), total_time=0.0)
        plt = plot_harv_mesh(case)
        @test plt isa Plots.Plot
    else
        @test_broken false
    end
end
