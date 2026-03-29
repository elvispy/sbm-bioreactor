using SBM_Bioreactor

function run_case(n; blocked::Bool, blocked_linear_solver::Bool=false, blocked_outer_solver=:gmres)
    case = build_harv_2d_case(
        partition = (n, n),
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

function warm_and_time(n; name, blocked, blocked_linear_solver=false, blocked_outer_solver=:gmres)
    first = run_case(n; blocked=blocked, blocked_linear_solver=blocked_linear_solver, blocked_outer_solver=blocked_outer_solver)
    second = run_case(n; blocked=blocked, blocked_linear_solver=blocked_linear_solver, blocked_outer_solver=blocked_outer_solver)
    println((name=name, n=n, first=first, second=second))
    flush(stdout)
end

function main()
    for n in (32, 64, 128)
        warm_and_time(n; name="blocked_gmres", blocked=true, blocked_linear_solver=true, blocked_outer_solver=:gmres)
        warm_and_time(n; name="unblocked_direct", blocked=false, blocked_linear_solver=false)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
