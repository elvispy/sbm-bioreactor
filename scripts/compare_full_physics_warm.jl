using SBM_Bioreactor

function run_case(; blocked::Bool, blocked_linear_solver::Bool=false, blocked_outer_solver=:gmres)
    case = build_harv_2d_case(
        partition = (12, 12),
        dt = 0.05,
        total_time = 0.05,
        degree = 2,
        blocked = blocked,
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
        blocked_linear_solver = blocked_linear_solver,
        blocked_outer_solver = blocked_outer_solver,
        transport_block_solver = :lu,
        profile_steps = true,
    )
    return result.profile.steps[1]
end

function warm_and_time(; name, blocked, blocked_linear_solver=false, blocked_outer_solver=:gmres)
    first = run_case(; blocked=blocked, blocked_linear_solver=blocked_linear_solver, blocked_outer_solver=blocked_outer_solver)
    second = run_case(; blocked=blocked, blocked_linear_solver=blocked_linear_solver, blocked_outer_solver=blocked_outer_solver)
    println((name=name, first=first, second=second))
end

function main()
    warm_and_time(; name="full_blocked_gmres", blocked=true, blocked_linear_solver=true, blocked_outer_solver=:gmres)
    warm_and_time(; name="full_blocked_direct", blocked=true, blocked_linear_solver=false)
    warm_and_time(; name="full_unblocked_direct", blocked=false, blocked_linear_solver=false)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
