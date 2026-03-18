using Gridap
using LineSearches: BackTracking

function coupled_bioreactor_residual(x, x_prevs, y, dt, params, order=1, t=0.0)
    u, p, Φ, C, Γ = x
    v, q, w, z, v_γ = y
    
    # Time derivative term (BDF1 or BDF2)
    if order == 1
        u_n, p_n, Φ_n, C_n, Γ_n = x_prevs[1]
        u_dot = (u - u_n) / dt
        Φ_dot = (Φ - Φ_n) / dt
        C_dot = (C - C_n) / dt
    else
        u_n, p_n, Φ_n, C_n, Γ_n = x_prevs[1]
        u_nn, p_nn, Φ_nn, C_nn, Γ_nn = x_prevs[2]
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
    L = params.L
    u_wall = params.u_wall
    kc = params.kc
    ke = params.ke
    d0 = params.d0
    
    # Density and viscosity (Krieger-Dougherty)
    ρ = (1.0 - Φ) * ρf + Φ * ρs
    visc_op(phi) = krieger_viscosity(phi; μf=μf, Φmax=Φmax)
    μ = visc_op ∘ Φ
    
    # Analytical derivative of μ with respect to Φ
    dμ_dΦ_op(phi) = 2.5 * μf * (1.0 - phi/Φmax)^(-2.5*Φmax - 1.0)
    ∇μ = (dμ_dΦ_op ∘ Φ) * ∇(Φ)
    
    # 5th variable: Shear Rate Projection (Smoothing)
    res_gamma = v_γ * (Γ - shear_rate(u))

    # Navier Stokes weak form with Hele-Shaw drag (B.6)
    drag_coeff = 4.0 * μf / (L^2)
    res_ns = (ρ * u_dot ⋅ v) + navier_stokes_weak_form(u, p, v, q, μ, ρ, g, u_wall, drag_coeff)
    
    # The continuity equation is modified
    flux = particle_flux(u, Φ, ∇(Φ), μ, ∇μ, a, ρs, ρf, μf, Φavg, g, Γ, ∇(Γ))
    res_continuity_rhs = ∇(q) ⋅ (flux * ((ρs - ρf) / (ρs * ρf)))
    
    # Cell Transport & Growth kinetics (B.7)
    source_phi = (π/6.0 * a^3) * kc * C * d0 * exp(ke * t)
    res_phi = (w * Φ_dot) + (w * (u ⋅ ∇(Φ))) + (∇(w) ⋅ (flux * ((ρs - ρf) / (ρs * ρf)))) - (w * source_phi)
    
    # Nutrient Transport
    rc = -kc * (Φ / (π/6.0 * a^3))
    res_C = (z * C_dot) + (z * (u ⋅ ∇(C))) + (Df * ∇(z) ⊙ ∇(C)) - (z * rc)
    
    return res_ns + res_continuity_rhs + res_phi + res_C + res_gamma
end

function run_bioreactor_simulation(X, Y, dΩ, dt, params, nsteps; write_vtk_interval=1)
    # Initial state
    x_n = interpolate_everywhere([params.u0, params.p0, params.Φ0, params.C0, params.Γ0], X)
    x_nn = x_n # For BDF2, first step uses BDF1
    
    # Newton Solver
    nls = NLSolver(show_trace=true, method=:newton, linesearch=BackTracking())
    solver = FESolver(nls)
    
    xh = x_n
    
    for step in 1:nsteps
        t = step * dt
        println("Step: $step, Time: $t")
        
        # Determine order
        order = step == 1 ? 1 : 2
        x_prevs = order == 1 ? (x_n,) : (x_n, x_nn)
        
        res(x, y) = ∫( coupled_bioreactor_residual(x, x_prevs, y, dt, params, order, t) )dΩ
        op = FEOperator(res, X, Y)
        
        xh, _ = solve!(xh, solver, op)
        
        # Update previous steps
        x_nn = x_n
        x_n = xh
        
        if step % write_vtk_interval == 0
            writevtk(get_triangulation(dΩ), "results_$step", cellfields=["u"=>xh[1], "p"=>xh[2], "phi"=>xh[3], "C"=>xh[4], "gamma"=>xh[5]])
        end
    end
    
    return xh
end
