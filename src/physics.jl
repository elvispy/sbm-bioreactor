function navier_stokes_weak_form(u, p, v, q, μ, ρ, f)
    # The convective term (ρ * (u ⋅ ∇(u)) ⋅ v), the viscous term (μ * ∇(u) ⊙ ∇(v)),
    # the pressure gradient (−p * (∇ ⋅ v)), the continuity equation (q * (∇ ⋅ u)),
    # and the body force term (−f ⋅ v).
    (ρ * (u ⋅ ∇(u)) ⋅ v) + (μ * ∇(u) ⊙ ∇(v)) - (p * (∇ ⋅ v)) + (q * (∇ ⋅ u)) - (f ⋅ v)
end

function krieger_viscosity(Φ; μf=0.5889, Φmax=0.64)
    return μf * (1 - Φ/Φmax)^(-2.5*Φmax)
end

