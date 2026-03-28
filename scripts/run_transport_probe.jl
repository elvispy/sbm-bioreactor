using SBM_Bioreactor
using Gridap
using Gridap.Algebra
using GridapSolvers
using GridapSolvers.LinearSolvers
using GridapSolvers.BlockSolvers
using LinearAlgebra

include("experimental_block_solver_linear_probe.jl")

function run_probe(nsizes::Vector{Int}, outfile::String)
    open(outfile, "w") do io
        for n in nsizes
            case, _, rhs, J = build_block_problem(n=n)
            println(io, (partition=case.metadata.partition, ndofs=num_free_dofs(case.X)))
            for transport_kind in (:lu, :gmres_jacobi)
                x, t = solve_block_iterative(J, rhs; transport_kind=transport_kind)
                rnorm, rhsnorm = residual_metrics(J, x, rhs)
                println(io, (
                    transport=String(transport_kind),
                    time=t,
                    residual_norm=rnorm,
                    relative_residual=rnorm / rhsnorm,
                    xnorm=norm(x),
                ))
                flush(io)
            end
        end
    end
end

function main()
    nsizes = isempty(ARGS) ? [64] : parse.(Int, ARGS[1:end-1])
    outfile = isempty(ARGS) ? "/tmp/transport_probe.txt" : ARGS[end]
    run_probe(nsizes, outfile)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
