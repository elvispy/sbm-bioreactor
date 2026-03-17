using Gridap

function coupled_bioreactor_residual(x, x_n, y, dt, params)
    u, p, Φ, C = x
    u_n, p_n, Φ_n, C_n = x_n
    v, q, w, z = y
    
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
    res_ns = (ρ * ((u - u_n) / dt) ⋅ v) + navier_stokes_weak_form(u, p, v, q, μ, ρ, g)
    
    # The continuity equation is modified
    flux = particle_flux(u, Φ, ∇(Φ), μ, ∇μ, a, ρs, ρf, μf, Φavg, g)
    # The pressure test function q tests the continuity: q * (∇⋅u - ... )
    # But navier_stokes_weak_form already has q * ∇⋅u, so we add the RHS:
    # Actually navier_stokes_weak_form returns: ... + q * ∇⋅u
    # We subtract the flux divergence: - q * ((ρs - ρf) / (ρs * ρf)) * ∇⋅flux
    # Or, integrating by parts: + ∇(q) ⋅ (flux * ((ρs - ρf) / (ρs * ρf)))
    res_continuity_rhs = ∇(q) ⋅ (flux * ((ρs - ρf) / (ρs * ρf)))
    
    # Cell Transport: ∂Φ/∂t + u⋅∇Φ = (ρs - ρf)/(ρs * ρf) ∇⋅Js
    # Integrated against test function w
    res_phi = (w * ((Φ - Φ_n) / dt)) + (w * (u ⋅ ∇(Φ))) + (∇(w) ⋅ (flux * ((ρs - ρf) / (ρs * ρf))))
    
    # Nutrient Transport: ∂C/∂t + u⋅∇C = Df ∇²C + rc
    # We ignore rc for the basic solver test
    res_C = (z * ((C - C_n) / dt)) + (z * (u ⋅ ∇(C))) + (Df * ∇(z) ⊙ ∇(C))
    
    return res_ns + res_continuity_rhs + res_phi + res_C
end

function run_bioreactor_step(model, degree, dt, params, u_n, p_n, Φ_n, C_n, U, P, Q_phi, Q_C)
    # This is a skeleton. A real implementation needs the FESpaces properly set up.
    pass
end
