"""
    navier_stokes_weak_form(u, p, v, q, μ, ρ, f, u_drag=nothing, drag_coeff=0.0; include_convection=true)

Compute the weak form of the incompressible Navier-Stokes equations with an optional 
Hele-Shaw depth-averaged friction term (B.6 in Chao & Das 2015).

# Arguments
- `u`, `p`: Trial velocity and pressure fields.
- `v`, `q`: Test velocity and pressure fields.
- `μ`: Dynamic viscosity (typically depends on particle concentration Φ).
- `ρ`: Mixture density.
- `f`: External body forces (e.g., gravity).
- `u_drag`: Reference velocity for drag (e.g., wall velocity).
- `drag_coeff`: Friction coefficient for depth-averaged flow.

# Physical Meaning
1. `(ρ * (u ⋅ ∇(u)) ⋅ v)`: Inertial term (convection of momentum).
2. `(μ * ∇(u) ⊙ ∇(v))`: Viscous stress term (dissipation).
3. `-(p * (∇ ⋅ v))`: Pressure gradient term.
4. `(q * (∇ ⋅ u))`: Continuity constraint (incompressibility).
5. `-(f ⋅ v)`: Body force term (e.g., buoyancy).
6. `(drag_coeff * (u - u_drag) ⋅ v)`: Hele-Shaw drag representing the out-of-plane 
   viscous resistance in thin-gap bioreactors.
"""
function navier_stokes_weak_form(u, p, v, q, μ, ρ, f, u_drag=nothing, drag_coeff=0.0; include_convection=true)
    # Standard Navier-Stokes terms: Advection + Diffusion - Pressure + Divergence Constraint - Forces
    convection = include_convection ? (ρ * (u ⋅ ∇(u)) ⋅ v) : (0.0 * (μ * ∇(u) ⊙ ∇(v)))
    res = convection + (μ * ∇(u) ⊙ ∇(v)) - (p * (∇ ⋅ v)) + (q * (∇ ⋅ u)) - (f ⋅ v)
    
    # Hele-Shaw Depth-Averaged Friction (B.6)
    # Models the viscous drag from the top and bottom plates in a 2D depth-averaged simulation.
    if u_drag !== nothing && drag_coeff > 0.0
        res = res + (drag_coeff * (u - u_drag) ⋅ v)
    end
    
    return res
end

"""
    krieger_viscosity(Φ; μf=0.5889, Φmax=0.64)

Calculate the effective mixture viscosity using the Krieger-Dougherty empirical relation.

# Arguments
- `Φ`: Local particle volume fraction.
- `μf`: Viscosity of the suspending fluid (base fluid).
- `Φmax`: Maximum packing fraction where viscosity diverges (typically ~0.64 for spheres).

# Physical Meaning
Viscosity increases nonlinearly with particle concentration, becoming infinite as Φ 
approaches Φmax, reflecting the transition from a fluid-like to a solid-like state.
"""
function krieger_viscosity(Φ; μf=0.5889, Φmax=0.64)
    return μf * (1 - Φ/Φmax)^(-2.5*Φmax)
end

"""
    shear_rate(u)

Compute the local shear rate (scalar invariant of the strain rate tensor).

# Arguments
- `u`: Velocity field.

# Physical Meaning
The shear rate Γ = sqrt(2 * ε(u) : ε(u)) quantifies the intensity of local fluid 
deformation, which drives the shear-induced migration of particles.
"""
function shear_rate(u)
    sqrt_op(x) = sqrt(abs(x) + 1e-10) # Regularized square root to avoid singularity at zero
    return sqrt_op ∘ (2.0 * ε(u) ⊙ ε(u))
end

"""
    particle_flux(u, Φ, ∇Φ, μ, ∇μ, a, ρs, ρf, μf, Φavg, g, Γ, ∇Γ)

Compute the total particle migration flux based on the Suspension Balance Model (SBM).

# Arguments
- `u`, `Φ`, `∇Φ`: Velocity, volume fraction, and its gradient.
- `μ`, `∇μ`: Mixture viscosity and its gradient.
- `a`: Particle radius.
- `ρs`, `ρf`: Particle and fluid densities.
- `μf`: Fluid viscosity.
- `Φavg`: Average volume fraction (used in sedimentation hindrance).
- `g`: Gravity vector.
- `Γ`, `∇Γ`: Local shear rate and its gradient.

# Physical Meaning (Chao & Das 2015)
1. `Jsc = -0.41 * a^2 * Φ * ∇(ΓΦ)`: Shear-induced migration from high to low shear 
   and concentration gradients (particle-particle collisions).
2. `Jsμ = -0.62 * a^2 * Φ^2 * Γ * ∇ln(μ)`: Migration toward lower viscosity regions.
3. `Jst = - (ust) * Φ`: Sedimentation flux due to buoyancy, corrected for hindrance 
   effects at high Φ.
"""
function particle_flux(u, Φ, ∇Φ, μ, ∇μ, a, ρs, ρf, μf, Φavg, g, Γ, ∇Γ)
    # Total gradient of (Γ * Φ) for the Jsc term: ∇(Γ*Φ) = Γ*∇Φ + Φ*∇Γ
    ∇ΓΦ = Γ * ∇Φ + Φ * ∇Γ

    # 1. Flux due to shear rate and concentration gradients (Chao & Das Eq. 12)
    Jsc = -0.41 * (a^2) * Φ * ∇ΓΦ

    if μ isa Number
        # Frozen viscosity removes the viscosity-gradient migration term entirely.
        ust = (μf * (1.0 - Φavg) / μ) * (2.0 * a^2 * (ρs - ρf) / (9.0 * μ)) * g
        Jst = -(ust * Φ)
        return Jsc + Jst
    end

    # Gradient of log-viscosity for the Jsμ term
    ∇lnμ = (1.0 / μ) * ∇μ

    # 2. Flux due to viscosity gradients (Chao & Das Eq. 13)
    Jsμ = -0.62 * (a^2) * (Φ * Φ) * Γ * ∇lnμ

    # 3. Sedimentation / Buoyancy flux (Chao & Das Eq. 14)
    # Corrected for fluid viscosity and hindrance via Richardson-Zaki like term
    ust_fh_op(m) = (μf * (1.0 - Φavg) / m) * (2.0 * a^2 * (ρs - ρf) / (9.0 * m)) * g
    ust = ust_fh_op ∘ μ
    Jst = -(ust * Φ)
    
    return Jsc + Jsμ + Jst
end
