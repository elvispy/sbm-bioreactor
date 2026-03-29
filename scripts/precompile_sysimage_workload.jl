using SBM_Bioreactor

function run_workload()
    case = build_harv_2d_case(
        partition = (8, 8),
        dt = 0.05,
        total_time = 0.10,
        degree = 2,
        blocked = true,
    )

    params = merge(case.params, (use_explicit_jacobian = true,))

    result = run_bioreactor_simulation(
        case.X,
        case.Y,
        case.dΩ,
        case.metadata.dt,
        params,
        1;
        write_vtk_interval = 0,
        nonlinear_show_trace = false,
        max_order = 1,
        blocked_linear_solver = true,
        blocked_outer_solver = :gmres,
        transport_block_solver = :lu,
        profile_steps = false,
    )

    @assert !isempty(result.profile.steps)
    return nothing
end

run_workload()
