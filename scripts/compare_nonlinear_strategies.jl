using SBM_Bioreactor
using Gridap
using Gridap.FESpaces: get_algebraic_operator
using SparseArrays
using LinearAlgebra
using LineSearches: BackTracking
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

function make_easy_case(; partition=(1, 1), dt=0.2, degree=2)
    case = build_harv_2d_case(partition=partition, dt=dt, total_time=dt, degree=degree)
    params = merge(
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
    x0_fe = interpolate_everywhere([params.u0, params.p0, params.Φ0, params.C0, params.Γ0], case.X)
    x_prevs = (x0_fe,)
    res(x, y) = ∫(coupled_bioreactor_residual(x, x_prevs, y, case.metadata.dt, params, 1, case.metadata.dt))case.dΩ
    jac(x, dx, y) = ∫(easy_manual_jacobian(x, x_prevs, dx, y, case.metadata.dt, params, 1, case.metadata.dt))case.dΩ
    op_fe = FEOperator(res, jac, case.X, case.Y)
    op = get_algebraic_operator(op_fe)
    x0 = copy(get_free_dof_values(x0_fe))
    block_lengths = map(num_free_dofs, getproperty(case.X, :spaces))
    return case, params, op, x0, block_lengths
end

function make_df(op, x0)
    f!(r, x) = Gridap.Algebra.residual!(r, op, x)
    j!(J, x) = Gridap.Algebra.jacobian!(J, op, x)
    fj!(r, J, x) = Gridap.Algebra.residual_and_jacobian!(r, J, op, x)
    f0, J0 = Gridap.Algebra.residual_and_jacobian(op, x0)
    return Gridap.Algebra.OnceDifferentiable(f!, j!, fj!, x0, f0, J0)
end

function build_scaling_vector(block_lengths, params)
    ω = params.u_wall(VectorValue(0.0, 1.0))[1]
    field_scales = [
        max(abs(ω), 1.0e-2),
        1.0,
        0.1,
        5.5,
        1.0e-5,
    ]
    parts = [fill(field_scales[i], block_lengths[i]) for i in eachindex(block_lengths)]
    return reduce(vcat, parts), field_scales
end

function make_scaled_df(op, x0, scales)
    z0 = x0 ./ scales
    r0, J0_sparse = Gridap.Algebra.residual_and_jacobian(op, x0)
    J0 = Matrix(J0_sparse) * Diagonal(scales)

    function f!(r, z)
        x = scales .* z
        tmp = similar(r0)
        Gridap.Algebra.residual!(tmp, op, x)
        copyto!(r, tmp)
    end

    function j!(J, z)
        x = scales .* z
        tmp = similar(J0_sparse)
        Gridap.Algebra.jacobian!(tmp, op, x)
        J[:, :] = Matrix(tmp) * Diagonal(scales)
    end

    function fj!(r, J, z)
        x = scales .* z
        tmp_r, tmp_J_sparse = Gridap.Algebra.residual_and_jacobian(op, x)
        copyto!(r, tmp_r)
        J[:, :] = Matrix(tmp_J_sparse) * Diagonal(scales)
    end

    return Gridap.Algebra.OnceDifferentiable(f!, j!, fj!, z0, r0, J0), z0
end

function print_result(name, result, elapsed; extra="")
    @printf("%s\n", name)
    @printf("  iterations = %d\n", result.iterations)
    @printf("  f_converged = %s\n", result.f_converged)
    @printf("  x_converged = %s\n", result.x_converged)
    @printf("  residual_norm = %.6e\n", result.residual_norm)
    @printf("  elapsed = %.3f s\n", elapsed)
    isempty(extra) || println(extra)
end

function run_strategy!(strategy, df, x0, op, scales; maxiters=200)
    println("running strategy=$(strategy), maxiters=$(maxiters)")
    flush(stdout)

    if strategy == "baseline"
        elapsed = @elapsed begin
            global result = Gridap.Algebra.nlsolve(
                df,
                copy(x0);
                method = :newton,
                linesearch = BackTracking(),
                show_trace = false,
                iterations = maxiters,
            )
        end
        print_result("baseline_newton_backtracking", result, elapsed, extra="  maxiters = $(maxiters)")
        return
    end

    if strategy == "xtol"
        elapsed = @elapsed begin
            global result = Gridap.Algebra.nlsolve(
                df,
                copy(x0);
                method = :newton,
                linesearch = BackTracking(),
                show_trace = false,
                iterations = maxiters,
                xtol = 1.0e-10,
            )
        end
        print_result("newton_with_xtol", result, elapsed, extra="  xtol = 1e-10, maxiters = $(maxiters)")
        return
    end

    if strategy == "trust_region"
        elapsed = @elapsed begin
            global result = Gridap.Algebra.nlsolve(
                df,
                copy(x0);
                method = :trust_region,
                show_trace = false,
                iterations = maxiters,
                autoscale = true,
            )
        end
        print_result("trust_region_autoscale", result, elapsed, extra="  autoscale = true, maxiters = $(maxiters)")
        return
    end

    if strategy == "scaled"
        scaled_df, z0 = make_scaled_df(op, x0, scales)
        elapsed = @elapsed begin
            global result = Gridap.Algebra.nlsolve(
                scaled_df,
                copy(z0);
                method = :newton,
                linesearch = BackTracking(),
                show_trace = false,
                iterations = maxiters,
            )
        end
        print_result("newton_with_field_scaling", result, elapsed, extra="  scaling acts on unknown blocks only, maxiters = $(maxiters)")
        return
    end

    error("Unknown strategy: $(strategy)")
end

partition = length(ARGS) >= 1 ? (parse(Int, ARGS[1]), parse(Int, ARGS[1])) : (1, 1)
degree = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2
strategy = length(ARGS) >= 3 ? ARGS[3] : "all"
maxiters = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 200

case, params, op, x0, block_lengths = make_easy_case(partition=partition, degree=degree)
scales, field_scales = build_scaling_vector(block_lengths, params)

println("nonlinear strategy comparison")
println("partition=$(partition), degree=$(degree), dofs=$(length(x0))")
println("block_lengths=$(block_lengths)")
println("field_scales=$(field_scales)")
println("requested_strategy=$(strategy)")
println("maxiters=$(maxiters)")

df = make_df(op, x0)

if strategy == "all"
    for name in ("baseline", "xtol", "trust_region", "scaled")
        run_strategy!(name, df, x0, op, scales; maxiters=maxiters)
    end
else
    run_strategy!(strategy, df, x0, op, scales; maxiters=maxiters)
end
