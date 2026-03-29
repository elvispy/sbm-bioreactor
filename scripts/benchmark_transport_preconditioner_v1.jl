using SBM_Bioreactor

function run_linear_probe(n; transport_kind::Symbol)
    case = build_harv_2d_case(
        partition = (n, n),
        dt = 0.05,
        total_time = 0.05,
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
        transport_block_solver = transport_kind,
        profile_steps = true,
    )
    return result.profile.steps[1]
end

function run_live_probe(; transport_kind::Symbol)
    case = build_harv_2d_case(
        partition = (12, 12),
        dt = 0.05,
        total_time = 0.05,
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
        transport_block_solver = transport_kind,
        profile_steps = true,
    )
    return result.profile.steps[1]
end

function warm_linear_probe(n; transport_kind::Symbol)
    first = run_linear_probe(n; transport_kind=transport_kind)
    second = run_linear_probe(n; transport_kind=transport_kind)
    return (first=first, second=second)
end

function main()
    for n in (64, 128)
        println((kind = :lu, n = n, result = warm_linear_probe(n; transport_kind=:lu)))
        flush(stdout)
    end
    println((kind = :lu, nonlinear_probe = run_live_probe(; transport_kind=:lu)))
    flush(stdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
