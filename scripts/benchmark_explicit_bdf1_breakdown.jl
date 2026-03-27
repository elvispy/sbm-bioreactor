using SBM_Bioreactor
using Gridap
using Gridap.Algebra: residual_and_jacobian
using SparseArrays

function build_case_data(n; degree=2, dt=0.1)
    case = build_harv_2d_case(partition=(n, n), dt=dt, total_time=dt, degree=degree)
    params = merge(case.params, (use_explicit_jacobian=true,))
    x0 = interpolate_everywhere([params.u0, params.p0, params.Φ0, params.C0, params.Γ0], case.X)
    x0_prevs = (x0,)
    op = build_bioreactor_operator(case.X, case.Y, case.dΩ, x0_prevs, case.metadata.dt, params, 1, case.metadata.dt)
    return case, params, x0, op
end

function warm_and_time_assembly(op, x0)
    residual_and_jacobian(op, x0)
    r = nothing
    J = nothing
    t = @elapsed begin
        r, J = residual_and_jacobian(op, x0)
    end
    return t, r, J
end

function warm_and_time_solve(case, params)
    run_bioreactor_simulation(
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
    timed = run_bioreactor_simulation(
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
    return timed.profile.steps[1].solve_time
end

function benchmark_case(n; degree=2, dt=0.1)
    println(("case_start", n))
    case, params, x0, op = build_case_data(n; degree=degree, dt=dt)
    ndofs = num_free_dofs(case.X)
    assembly_time, r, J = warm_and_time_assembly(op, x0)
    solve_time = warm_and_time_solve(case, params)
    println((
        n = n,
        ndofs = ndofs,
        residual_length = length(r),
        nnz = nnz(J),
        assembly_time = assembly_time,
        solve_time = solve_time,
    ))
end

for n in (64, 128, 256)
    benchmark_case(n)
end
