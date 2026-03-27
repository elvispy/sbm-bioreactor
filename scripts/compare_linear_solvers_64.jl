using SBM_Bioreactor
using Gridap
using GridapPETSc
using LineSearches

function build_problem(; n=64, degree=2, dt=0.1)
    case = build_harv_2d_case(partition=(n, n), dt=dt, total_time=dt, degree=degree)
    params = merge(case.params, (use_explicit_jacobian=true,))
    x0 = interpolate_everywhere([params.u0, params.p0, params.Φ0, params.C0, params.Γ0], case.X)
    x0_prevs = (x0,)
    op = SBM_Bioreactor.build_bioreactor_operator(case.X, case.Y, case.dΩ, x0_prevs, case.metadata.dt, params, 1, case.metadata.dt)
    return case, params, x0, op
end

function timed_one_step(ls, case, x0, op; label)
    nls = NLSolver(
        ls;
        show_trace=false,
        method=:newton,
        iterations=1000,
        xtol=0.0,
        ftol=1.0e-8,
        linesearch=BackTracking(),
    )
    solver = FESolver(nls)

    warm_guess = deepcopy(x0)
    solve!(warm_guess, solver, op)

    timed_guess = deepcopy(x0)
    t = @elapsed solve!(timed_guess, solver, op)
    println((label=label, solve_time=t))
end

case, params, x0, op = build_problem()
println((ndofs=num_free_dofs(case.X), partition=case.metadata.partition))

timed_one_step(BackslashSolver(), case, x0, op; label="direct_backslash")

options = "-ksp_type gmres -ksp_rtol 1.0e-8 -ksp_atol 1.0e-12 -pc_type ilu -ksp_converged_reason -ksp_error_if_not_converged true"
GridapPETSc.with(args=split(options)) do
    timed_one_step(PETScLinearSolver(), case, x0, op; label="petsc_gmres_ilu")
end
