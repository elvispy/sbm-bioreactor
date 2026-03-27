using SBM_Bioreactor

function run_case(n; degree=2, dt=0.2)
    case = build_harv_2d_case(partition=(n, n), dt=dt, total_time=dt, degree=degree)
    params = merge(case.params, (use_explicit_jacobian=true,))
    return run_bioreactor_simulation(
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
end

function benchmark_case(n; degree=2, dt=0.2)
    println(("warmup_start", n))
    warm = run_case(n; degree=degree, dt=dt)
    warm_step = warm.profile.steps[1]
    println((
        warmup = n,
        initial_setup_time = warm.profile.initial_setup_time,
        operator_build_time = warm_step.operator_build_time,
        solve_time = warm_step.solve_time,
    ))

    println(("timed_start", n))
    timed = run_case(n; degree=degree, dt=dt)
    step = timed.profile.steps[1]
    println((
        timed = n,
        initial_setup_time = timed.profile.initial_setup_time,
        operator_build_time = step.operator_build_time,
        solve_time = step.solve_time,
    ))
end

for n in (4, 8, 12, 16)
    benchmark_case(n)
end
