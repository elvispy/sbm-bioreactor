function navier_stokes_weak_form(u, p, v, q, μ, ρ, f, u_drag=nothing, drag_coeff=0.0)
    # Standard Navier-Stokes terms
    res = (ρ * (u ⋅ ∇(u)) ⋅ v) + (μ * ∇(u) ⊙ ∇(v)) - (p * (∇ ⋅ v)) + (q * (∇ ⋅ u)) - (f ⋅ v)
    
    # Hele-Shaw Depth-Averaged Friction (B.6)
    if u_drag !== nothing && drag_coeff > 0.0
        res = res + (drag_coeff * (u - u_drag) ⋅ v)
    end
    
    return res
end

function krieger_viscosity(Φ; μf=0.5889, Φmax=0.64)
    return μf * (1 - Φ/Φmax)^(-2.5*Φmax)
end

function shear_rate(u)
    sqrt_op(x) = sqrt(abs(x) + 1e-10)
    return sqrt_op ∘ (2.0 * ε(u) ⊙ ε(u))
end

function particle_flux(u, Φ, ∇Φ, μ, ∇μ, a, ρs, ρf, μf, Φavg, g)
    γ̇ = shear_rate(u)
    
    inv_op(x) = 1.0 / x
    ∇lnμ = (inv_op ∘ μ) * ∇μ
    
    # We use ∇(γ̇*Φ) ≈ γ̇*∇Φ. 
    # Computing ∇(γ̇) requires 2nd derivatives of u which are not natively supported 
    # without a mixed formulation or projection to a continuous space. 
    # TODO: Implement L2 projection of γ̇ to a P1 space to compute ∇(γ̇).
    ∇γ̇Φ = γ̇ * ∇Φ

    Jsc = -0.41 * (a^2) * Φ * ∇γ̇Φ
    Jsμ = -0.62 * (a^2) * (Φ * Φ) * γ̇ * ∇lnμ
    
    ust_fh_op(m) = (μf * (1.0 - Φavg) / m) * (2.0 * a^2 * (ρs - ρf) / (9.0 * m)) * g
    Jst = - (ust_fh_op ∘ μ) * Φ
    
    return Jsc + Jsμ + Jst
end


