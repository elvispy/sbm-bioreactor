using SBM_Bioreactor
using Gridap
using Gridap.FESpaces: get_algebraic_operator
using SparseArrays
using LinearAlgebra
using LineSearches: BackTracking
using Printf

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
    op_fe = SBM_Bioreactor.build_bioreactor_operator(
        case.X,
        case.Y,
        case.dΩ,
        (x0_fe,),
        case.metadata.dt,
        params,
        1,
        case.metadata.dt,
    )
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
    @printf("  residual_inf = %.6e\n", maximum(abs, result.f))
    @printf("  elapsed = %.3f s\n", elapsed)
    isempty(extra) || println(extra)
end

partition = length(ARGS) >= 1 ? (parse(Int, ARGS[1]), parse(Int, ARGS[1])) : (1, 1)
degree = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2

case, params, op, x0, block_lengths = make_easy_case(partition=partition, degree=degree)
scales, field_scales = build_scaling_vector(block_lengths, params)

println("nonlinear strategy comparison")
println("partition=$(partition), degree=$(degree), dofs=$(length(x0))")
println("block_lengths=$(block_lengths)")
println("field_scales=$(field_scales)")

df = make_df(op, x0)

t_baseline = @elapsed begin
    global baseline = Gridap.Algebra.nlsolve(
        df,
        copy(x0);
        method = :newton,
        linesearch = BackTracking(),
        show_trace = false,
        iterations = 1000,
    )
end
print_result("baseline_newton_backtracking", baseline, t_baseline)

t_xtol = @elapsed begin
    global xtol_run = Gridap.Algebra.nlsolve(
        df,
        copy(x0);
        method = :newton,
        linesearch = BackTracking(),
        show_trace = false,
        iterations = 1000,
        xtol = 1.0e-10,
    )
end
print_result("newton_with_xtol", xtol_run, t_xtol, extra="  xtol = 1e-10")

t_tr = @elapsed begin
    global trust_region_run = Gridap.Algebra.nlsolve(
        df,
        copy(x0);
        method = :trust_region,
        show_trace = false,
        iterations = 1000,
        autoscale = true,
    )
end
print_result("trust_region_autoscale", trust_region_run, t_tr)

scaled_df, z0 = make_scaled_df(op, x0, scales)
t_scaled = @elapsed begin
    global scaled_run = Gridap.Algebra.nlsolve(
        scaled_df,
        copy(z0);
        method = :newton,
        linesearch = BackTracking(),
        show_trace = false,
        iterations = 1000,
    )
end
print_result("newton_with_field_scaling", scaled_run, t_scaled, extra="  scaling acts on unknown blocks only")
