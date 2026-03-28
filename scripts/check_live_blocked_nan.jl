using SBM_Bioreactor
using Gridap
using Gridap.Algebra
using GridapSolvers
using GridapSolvers.LinearSolvers
using GridapSolvers.BlockSolvers
using LinearAlgebra

function flatten_values(x)
    if x isa Number
        return Float64[x]
    elseif x isa AbstractArray
        vals = Float64[]
        for xi in x
            append!(vals, flatten_values(xi))
        end
        return vals
    else
        try
            return flatten_values(collect(x))
        catch
            return Float64[]
        end
    end
end

function stats(x)
    vals = flatten_values(x)
    return (
        count = length(vals),
        any_nan = any(isnan, vals),
        any_inf = any(isinf, vals),
        maxabs = isempty(vals) ? 0.0 : maximum(abs, vals),
    )
end

function build_block_solver()
    flow_solver = LUSolver()
    transport_solver = LUSolver()
    blocks = [
        LinearSystemBlock() LinearSystemBlock()
        LinearSystemBlock() LinearSystemBlock()
    ]
    coeffs = [
        1.0 1.0
        0.0 1.0
    ]
    preconditioner = BlockTriangularSolver(blocks, [flow_solver, transport_solver], coeffs, :upper)
    return FGMRESSolver(
        20,
        preconditioner;
        maxiter = 100,
        atol = 1.0e-12,
        rtol = 1.0e-8,
        verbose = false,
    )
end

function main()
    case = build_harv_2d_case(
        partition = (12, 12),
        dt = 0.05,
        total_time = 0.05,
        degree = 2,
        blocked = true,
    )
    params = merge(
        case.params,
        (
            use_explicit_jacobian = true,
            freeze_viscosity = true,
            enable_particle_flux = false,
            enable_growth_source = false,
            enable_nutrient_reaction = false,
            include_convection = false,
        ),
    )
    x0 = interpolate_everywhere(
        [params.u0, params.p0, params.Φ0, params.C0, params.Γ0],
        case.X,
    )
    xfree = get_free_dof_values(x0)
    op = SBM_Bioreactor.build_bioreactor_operator(
        case.X,
        case.Y,
        case.dΩ,
        (x0,),
        case.metadata.dt,
        params,
        1,
        case.metadata.dt,
    )
    r, J = residual_and_jacobian(op, x0)
    println((rhs_type = typeof(r), matrix_type = typeof(J)))
    println((rhs_stats = stats(r), J_stats = stats(J)))
    for i in 1:2, j in 1:2
        blk = getfield(J, :blocks)[i, j]
        println((index = (i, j), type = typeof(blk), stats = stats(blk)))
    end
    for i in 1:2
        b = getfield(r, :blocks)[i]
        println((rhs_block = i, type = typeof(b), stats = stats(b)))
    end

    solver = build_block_solver()

    ss_plain = symbolic_setup(solver, J)
    ns_plain = numerical_setup(ss_plain, J)
    y_plain = allocate_in_domain(J)
    fill!(y_plain, 0.0)
    solve!(y_plain, ns_plain, r)
    println((plain_setup_solution_stats = stats(y_plain),))
    Jy_plain = allocate_in_range(J)
    fill!(Jy_plain, 0.0)
    mul!(Jy_plain, J, y_plain)
    println((plain_setup_norms = (rhs = norm(r), y = norm(y_plain), Jy = norm(Jy_plain)),))
    println((plain_setup_matvec_stats = stats(Jy_plain),))

    ss_x = symbolic_setup(solver, J, xfree)
    ns_x = numerical_setup(ss_x, J, xfree)
    y_x = allocate_in_domain(J)
    fill!(y_x, 0.0)
    solve!(y_x, ns_x, r)
    println((x_setup_solution_stats = stats(y_x),))
    Jy_x = allocate_in_range(J)
    fill!(Jy_x, 0.0)
    mul!(Jy_x, J, y_x)
    println((x_setup_norms = (rhs = norm(r), y = norm(y_x), Jy = norm(Jy_x)),))
    println((x_setup_matvec_stats = stats(Jy_x),))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
