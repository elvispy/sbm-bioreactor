using Gridap
using SBM_Bioreactor

function easy_manual_jacobian_guard(x, x_prevs, dx, y, dt, params, order=1, t=0.0)
    order == 1 || error("guard helper only supports BDF1")

    u, p, Φ, C, Γ = x
    du, dp, dΦ, dC, dΓ = dx
    v, q, w, z, v_γ = y
    u_n, _, _, _, _ = x_prevs[1]

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

case = build_harv_2d_case(partition=(1, 1), dt=0.2, total_time=0.2, degree=2)
x0 = interpolate_everywhere([case.params.u0, case.params.p0, case.params.Φ0, case.params.C0, case.params.Γ0], case.X)
x0_prevs = (x0,)

easy_params = merge(
    case.params,
    (
        enable_particle_flux = false,
        freeze_viscosity = true,
        include_convection = false,
        enable_growth_source = false,
        enable_nutrient_reaction = false,
        use_explicit_jacobian = true,
    ),
)

jac_easy_pkg = coupled_bioreactor_jacobian(x0, x0_prevs, x0, x0, case.metadata.dt, easy_params, 1, case.metadata.dt)
jac_easy_manual = easy_manual_jacobian_guard(x0, x0_prevs, x0, x0, case.metadata.dt, easy_params, 1, case.metadata.dt)
easy_diff = sum(∫(jac_easy_pkg - jac_easy_manual)case.dΩ)
println(("easy_diff", easy_diff))
abs(easy_diff) <= 1.0e-10 || error("easy explicit Jacobian no longer matches manual reference")

richer_params = merge(
    case.params,
    (
        enable_particle_flux = false,
        freeze_viscosity = true,
        include_convection = true,
        enable_growth_source = true,
        enable_nutrient_reaction = true,
        use_explicit_jacobian = true,
    ),
)

richer_op = SBM_Bioreactor.build_bioreactor_operator(
    case.X,
    case.Y,
    case.dΩ,
    x0_prevs,
    case.metadata.dt,
    richer_params,
    1,
    case.metadata.dt,
)
println(("richer_operator", typeof(richer_op)))
