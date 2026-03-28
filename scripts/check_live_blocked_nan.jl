using SBM_Bioreactor
using Gridap
using Gridap.Algebra

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
        blk = J.array[i, j]
        println((index = (i, j), type = typeof(blk), stats = stats(blk)))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
