using SBM_Bioreactor

function benchmark_case(n; degree=2, dt=0.2)
    case = build_harv_2d_case(partition=(n, n), dt=dt, total_time=dt, degree=degree)
    params = merge(case.params, (use_explicit_jacobian=true,))
    result = run_bioreactor_simulation(
        case.X,
        case.Y,
        case.dΩ,
        case.metadata.dt,
        params,
        1;
        write_vtk_interval=0,
        profile_steps=true,
        nonlinear_show_trace=false,
    )
    step = result.profile.steps[1]
    return (
        n = n,
        initial_setup_time = result.profile.initial_setup_time,
        operator_build_time = step.operator_build_time,
        solve_time = step.solve_time,
    )
end

for n in (1, 2, 4)
    println(("benchmark_start", n))
    stats = benchmark_case(n)
    println(stats)
end
