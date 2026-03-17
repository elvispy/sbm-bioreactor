using Gridap
using LineSearches: BackTracking

function coupled_bioreactor_residual(x, x_prevs, y, dt, params, order=1)
    u, p, Φ, C = x
    v, q, w, z = y
    
    # Time derivative term (BDF1 or BDF2)
    if order == 1
        # BDF1 (Backward Euler)
        u_n, p_n, Φ_n, C_n = x_prevs[1]
        u_dot = (u - u_n) / dt
        Φ_dot = (Φ - Φ_n) / dt
        C_dot = (C - C_n) / dt
    else
        # BDF2
        u_n, p_n, Φ_n, C_n = x_prevs[1]
        u_nn, p_nn, Φ_nn, C_nn = x_prevs[2]
        u_dot = (3.0*u - 4.0*u_n + u_nn) / (2.0*dt)
        Φ_dot = (3.0*Φ - 4.0*Φ_n + Φ_nn) / (2.0*dt)
        C_dot = (3.0*C - 4.0*C_n + C_nn) / (2.0*dt)
    end

    # Unpack params
    μf = params.μf
    Φmax = params.Φmax
    a = params.a
    ρs = params.ρs
    ρf = params.ρf
    g = params.g
    Df = params.Df
    Φavg = params.Φavg
    
    # Density and viscosity (Krieger-Dougherty)
    ρ = (1.0 - Φ) * ρf + Φ * ρs
    
    visc_op(phi) = krieger_viscosity(phi; μf=μf, Φmax=Φmax)
    μ = visc_op ∘ Φ
    
    # Analytical derivative of μ with respect to Φ
    dμ_dΦ_op(phi) = 2.5 * μf * (1.0 - phi/Φmax)^(-2.5*Φmax - 1.0)
    ∇μ = (dμ_dΦ_op ∘ Φ) * ∇(Φ)
    
    # Navier Stokes weak form
    res_ns = (ρ * u_dot ⋅ v) + navier_stokes_weak_form(u, p, v, q, μ, ρ, g)
    
    # The continuity equation is modified
    flux = particle_flux(u, Φ, ∇(Φ), μ, ∇μ, a, ρs, ρf, μf, Φavg, g)
    res_continuity_rhs = ∇(q) ⋅ (flux * ((ρs - ρf) / (ρs * ρf)))
    
    # Cell Transport: ∂Φ/∂t + u⋅∇Φ = (ρs - ρf)/(ρs * ρf) ∇⋅Js
    res_phi = (w * Φ_dot) + (w * (u ⋅ ∇(Φ))) + (∇(w) ⋅ (flux * ((ρs - ρf) / (ρs * ρf))))
    
    # Nutrient Transport: ∂C/∂t + u⋅∇C = Df ∇²C + rc
    res_C = (z * C_dot) + (z * (u ⋅ ∇(C))) + (Df * ∇(z) ⊙ ∇(C))
    
    return res_ns + res_continuity_rhs + res_phi + res_C
end

function run_bioreactor_simulation(X, Y, dΩ, dt, params, nsteps; write_vtk_interval=1)
    # Initial state
    x_n = interpolate_everywhere([params.u0, params.p0, params.Φ0, params.C0], X)
    x_nn = x_n # For BDF2, first step uses BDF1
    
    # Newton Solver
    nls = NLSolver(show_trace=true, method=:newton, linesearch=BackTracking())
    solver = FESolver(nls)
    
    xh = x_n
    
    for step in 1:nsteps
        println("Step: $step")
        
        # Determine order
        order = step == 1 ? 1 : 2
        x_prevs = order == 1 ? (x_n,) : (x_n, x_nn)
        
        res(x, y) = ∫( coupled_bioreactor_residual(x, x_prevs, y, dt, params, order) )dΩ
        op = FEOperator(res, X, Y)
        
        xh, _ = solve!(xh, solver, op)
        
        # Update previous steps
        x_nn = x_n
        x_n = xh
        
        if step % write_vtk_interval == 0
            writevtk(get_triangulation(dΩ), "results_$step", cellfields=["u"=>xh[1], "p"=>xh[2], "phi"=>xh[3], "C"=>xh[4]])
        end
    end
    
    return xh
end
