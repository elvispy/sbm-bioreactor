using SBM_Bioreactor
using Gridap
using Gridap.Algebra
using GridapSolvers
using GridapSolvers.LinearSolvers
using GridapSolvers.BlockSolvers
using LinearAlgebra

function build_block_problem(; n=64, degree=2, dt=0.1)
    case = build_harv_2d_case(
        partition=(n, n),
        dt=dt,
        total_time=dt,
        degree=degree,
        blocked=true,
    )
    params = merge(case.params, (use_explicit_jacobian=true,))
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
    return case, x0, rhs, J
end


function build_transport_solver(kind::Symbol)
    if kind == :gmres_jacobi
        return GMRESSolver(
            20;
            Pr=JacobiLinearSolver(),
            maxiter=200,
            atol=1.0e-12,
            rtol=1.0e-8,
            verbose=false,
        )
    elseif kind == :gmres_symgs
        return GMRESSolver(
            20;
            Pr=LinearSolverFromSmoother(SymGaussSeidelSmoother(2, 1.0)),
            maxiter=200,
            atol=1.0e-12,
            rtol=1.0e-8,
            verbose=false,
        )
    elseif kind == :gmres_richardson_symgs
        return GMRESSolver(
            20;
            Pr=RichardsonSmoother(SymGaussSeidelSmoother(1, 1.0), 3, 0.8),
            maxiter=200,
            atol=1.0e-12,
            rtol=1.0e-8,
            verbose=false,
        )
    else
        error("unsupported transport solver kind: $kind")
    end
end

function build_block_preconditioner(J; transport_kind::Symbol=:gmres_jacobi)
    flow_solver = LUSolver()
    transport_solver = build_transport_solver(transport_kind)
    blocks = [
        LinearSystemBlock() LinearSystemBlock()
        LinearSystemBlock() LinearSystemBlock()
    ]
    coeffs = [
        1.0 1.0
        0.0 1.0
    ]
    preconditioner = BlockTriangularSolver(blocks, [flow_solver, transport_solver], coeffs, :upper)
    solver = FGMRESSolver(
        20,
        preconditioner;
        maxiter=100,
        atol=1.0e-12,
        rtol=1.0e-8,
        verbose=true,
    )
    ns = numerical_setup(symbolic_setup(solver, J), J)
    return solver, ns
end

function solve_block_iterative(J, rhs; transport_kind::Symbol=:gmres_jacobi)
    _, ns = build_block_preconditioner(J; transport_kind=transport_kind)
    x = allocate_in_domain(J)
    fill!(x, 0.0)
    t = @elapsed solve!(x, ns, rhs)
    return x, t
end

function residual_metrics(J, x, rhs)
    Ax = allocate_in_range(J)
    mul!(Ax, J, x)
    r = similar(rhs)
    copy!(r, rhs)
    r .-= Ax
    return norm(r), norm(rhs)
end

function main(; n=64, degree=2, dt=0.1)
    case, _, rhs, J = build_block_problem(; n=n, degree=degree, dt=dt)
    println((partition=case.metadata.partition, degree=case.metadata.degree, blocked=case.metadata.blocked, ndofs=num_free_dofs(case.X)))
    println((matrix_type=string(typeof(J)), rhs_type=string(typeof(rhs))))

    for transport_kind in (:gmres_jacobi, :gmres_symgs, :gmres_richardson_symgs)
        try
            x_block, t_block = solve_block_iterative(J, rhs; transport_kind=transport_kind)
            rnorm, rhsnorm = residual_metrics(J, x_block, rhs)
            println((
                solver="block_fgmres_upper",
                transport=String(transport_kind),
                time=t_block,
                xnorm=norm(x_block),
                residual_norm=rnorm,
                relative_residual=rnorm / rhsnorm,
            ))
        catch err
            println((solver="block_fgmres_upper", transport=String(transport_kind), error=sprint(showerror, err)))
        end
    end
end

main()
