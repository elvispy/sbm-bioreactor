using SBM_Bioreactor
using Gridap
using Gridap.Algebra
using GridapSolvers
using GridapSolvers.LinearSolvers
using GridapSolvers.BlockSolvers
using LinearAlgebra

function build_case_and_system()
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
    rhs = similar(r)
    copy!(rhs, r)
    rmul!(rhs, -1)
    return J, rhs
end

function build_preconditioner()
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
    return BlockTriangularSolver(blocks, [flow_solver, transport_solver], coeffs, :upper)
end

function residual_metrics(J, x, rhs)
    Ax = allocate_in_range(J)
    fill!(Ax, 0.0)
    mul!(Ax, J, x)
    r = similar(rhs)
    copy!(r, rhs)
    r .-= Ax
    return norm(r), norm(rhs)
end

function run_solver(name, solver, J, rhs)
    println(("solver_start", name))
    ss = symbolic_setup(solver, J)
    ns = numerical_setup(ss, J)
    x = allocate_in_domain(J)
    fill!(x, 0.0)
    t = @elapsed solve!(x, ns, rhs)
    rnorm, rhsnorm = residual_metrics(J, x, rhs)
    println((
        solver = name,
        time = t,
        xnorm = norm(x),
        residual_norm = rnorm,
        relative_residual = rnorm / rhsnorm,
    ))
end

function main()
    J, rhs = build_case_and_system()
    P = build_preconditioner()
    gmres = GMRESSolver(20; Pr = P, maxiter = 100, atol = 1.0e-12, rtol = 1.0e-8, verbose = true)
    fgmres = FGMRESSolver(20, P; maxiter = 100, atol = 1.0e-12, rtol = 1.0e-8, verbose = true)
    run_solver("gmres_pr", gmres, J, rhs)
    run_solver("fgmres", fgmres, J, rhs)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
