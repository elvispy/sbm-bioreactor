using SBM_Bioreactor
using Gridap
using Gridap.FESpaces: residual_and_jacobian
using LineSearches: BackTracking
using LinearAlgebra
using Printf

function easy_manual_jacobian(x, x_prevs, dx, y, dt, params, order=1, t=0.0)
    order == 1 || error("Explicit easy-case Jacobian is implemented only for BDF1 in this experiment")

    u, p, Φ, C, Γ = x
    du, dp, dΦ, dC, dΓ = dx
    v, q, w, z, v_γ = y
    u_n, _, Φ_n, C_n, _ = x_prevs[1]

    μf = params.μf
    ρs = params.ρs
    ρf = params.ρf
    Df = params.Df
    L = params.L

    ρ = (1.0 - Φ) * ρf + Φ * ρs
    dρ = (ρs - ρf) * dΦ
    drag_coeff = 4.0 * μf / (L^2)

    u_dot = (u - u_n) / dt
    du_dot = du / dt
    dΦ_dot = dΦ / dt
    dC_dot = dC / dt

    shear_arg = 2.0 * ε(u) ⊙ ε(u)
    dshear_arg = 4.0 * ε(u) ⊙ ε(du)
    d_shear_rate_op(x) = sign(x) / (2.0 * sqrt(abs(x) + 1.0e-10))
    dshear = (d_shear_rate_op ∘ shear_arg) * dshear_arg

    jac_ns =
        (dρ * (u_dot ⋅ v)) +
        (ρ * (du_dot ⋅ v)) +
        (μf * ∇(du) ⊙ ∇(v)) -
        (dp * (∇ ⋅ v)) +
        (q * (∇ ⋅ du)) +
        (drag_coeff * (du ⋅ v))

    jac_phi =
        (w * dΦ_dot) +
        (w * (du ⋅ ∇(Φ))) +
        (w * (u ⋅ ∇(dΦ)))

    jac_C =
        (z * dC_dot) +
        (z * (du ⋅ ∇(C))) +
        (z * (u ⋅ ∇(dC))) +
        (Df * ∇(z) ⊙ ∇(dC))

    jac_gamma = v_γ * (dΓ - dshear)

    return jac_ns + jac_phi + jac_C + jac_gamma
end

function copy_state(X, xh)
    FEFunction(X, copy(get_free_dof_values(xh)))
end

partition = length(ARGS) >= 1 ? (parse(Int, ARGS[1]), parse(Int, ARGS[1])) : (2, 2)
degree = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2
mode = length(ARGS) >= 3 ? ARGS[3] : "both"
dt = 0.2

case = build_harv_2d_case(partition=partition, dt=dt, total_time=dt, degree=degree)
params = merge(
    case.params,
    (
        enable_particle_flux = false,
        freeze_viscosity = true,
        include_convection = false,
        enable_growth_source = false,
        enable_nutrient_reaction = false,
    ),
)
x0 = interpolate_everywhere([params.u0, params.p0, params.Φ0, params.C0, params.Γ0], case.X)
x_prevs = (x0,)

res(x, y) = ∫(coupled_bioreactor_residual(x, x_prevs, y, dt, params, 1, dt))case.dΩ
jac_explicit(x, dx, y) = ∫(easy_manual_jacobian(x, x_prevs, dx, y, dt, params, 1, dt))case.dΩ

op_ad = FEOperator(res, case.X, case.Y)
op_ex = FEOperator(res, jac_explicit, case.X, case.Y)

run_ad = mode in ("both", "autodiff")
run_ex = mode in ("both", "explicit")
run_ad || run_ex || error("mode must be one of: both, autodiff, explicit")

println("easy explicit-jacobian experiment")
println("partition=$(partition), degree=$(degree), dofs=$(length(get_free_dof_values(x0))), mode=$(mode)")

println("\nWarm-up")
if run_ad
    println("warming autodiff residual_and_jacobian ...")
    residual_and_jacobian(op_ad, x0)
end
if run_ex
    println("warming explicit residual_and_jacobian ...")
    residual_and_jacobian(op_ex, x0)
end
println("warm-up complete")

println("\nResidual/Jacobian agreement")
if run_ad && run_ex
    println("assembling autodiff residual/jacobian ...")
    b_ad, A_ad = residual_and_jacobian(op_ad, x0)
    println("assembling explicit residual/jacobian ...")
    b_ex, A_ex = residual_and_jacobian(op_ex, x0)
    db = b_ad - b_ex
    dA = Matrix(A_ad - A_ex)
    @printf("maxabs_residual_diff = %.6e\n", maximum(abs, db))
    @printf("l2_residual_diff = %.6e\n", norm(db))
    @printf("maxabs_jacobian_diff = %.6e\n", maximum(abs, dA))
    @printf("frobenius_jacobian_diff = %.6e\n", norm(dA))
else
    println("skipping agreement check because only one operator variant was requested")
end

println("\nFused residual_and_jacobian timing")
t_ad = NaN
t_ex = NaN
if run_ad
    println("timing autodiff residual_and_jacobian ...")
    t_ad = @elapsed residual_and_jacobian(op_ad, x0)
    @printf("autodiff_time = %.3f s\n", t_ad)
end
if run_ex
    println("timing explicit residual_and_jacobian ...")
    t_ex = @elapsed residual_and_jacobian(op_ex, x0)
    @printf("explicit_time = %.3f s\n", t_ex)
end
if run_ad && run_ex
    @printf("speedup = %.2fx\n", t_ad / t_ex)
end

println("\nOne Newton iteration timing")
nls = NLSolver(show_trace=true, method=:newton, iterations=1, linesearch=BackTracking())
solver = FESolver(nls)

t_one_ad = NaN
t_one_ex = NaN
t_one_ad_warm = NaN
t_one_ex_warm = NaN
if run_ad
    x_ad = copy_state(case.X, x0)
    t_one_ad = @elapsed begin
        println("running autodiff one-step solve ...")
        try
            solve!(x_ad, solver, op_ad)
        catch err
            println("autodiff_one_step_error=$(typeof(err))")
            println(err)
        end
    end
    @printf("autodiff_one_step_time = %.3f s\n", t_one_ad)
    x_ad_warm = copy_state(case.X, x0)
    t_one_ad_warm = @elapsed begin
        println("running autodiff warm one-step solve ...")
        try
            solve!(x_ad_warm, solver, op_ad)
        catch err
            println("autodiff_warm_one_step_error=$(typeof(err))")
            println(err)
        end
    end
    @printf("autodiff_warm_one_step_time = %.3f s\n", t_one_ad_warm)
end
if run_ex
    x_ex = copy_state(case.X, x0)
    t_one_ex = @elapsed begin
        println("running explicit one-step solve ...")
        try
            solve!(x_ex, solver, op_ex)
        catch err
            println("explicit_one_step_error=$(typeof(err))")
            println(err)
        end
    end
    @printf("explicit_one_step_time = %.3f s\n", t_one_ex)
    x_ex_warm = copy_state(case.X, x0)
    t_one_ex_warm = @elapsed begin
        println("running explicit warm one-step solve ...")
        try
            solve!(x_ex_warm, solver, op_ex)
        catch err
            println("explicit_warm_one_step_error=$(typeof(err))")
            println(err)
        end
    end
    @printf("explicit_warm_one_step_time = %.3f s\n", t_one_ex_warm)
end
if run_ad && run_ex
    @printf("one_step_speedup = %.2fx\n", t_one_ad / t_one_ex)
    @printf("warm_one_step_speedup = %.2fx\n", t_one_ad_warm / t_one_ex_warm)
end
