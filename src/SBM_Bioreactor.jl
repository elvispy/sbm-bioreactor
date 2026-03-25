"""
    module SBM_Bioreactor

A Julia package for simulating particle migration and cell growth in bioreactors 
using the Suspension Balance Model (SBM).

The solver implements a monolithic 5-variable finite element formulation (u, p, Φ, C, Γ) 
based on Gridap.jl. It accounts for shear-induced migration, buoyancy-driven 
sedimentation, and Hele-Shaw depth-averaged viscous effects as described in 
Chao & Das (2015).

# Key Features
- Monolithic coupling of Navier-Stokes and particle/nutrient transport.
- Krieger-Dougherty constitutive model for concentration-dependent viscosity.
- Shear-induced migration fluxes (SBM).
- Support for BDF1/BDF2 time stepping.
"""
module SBM_Bioreactor
using Gridap

include("physics.jl")
include("solver.jl")

"""run_bioreactor_simulation(X, Y, dΩ, dt, params, nsteps; write_vtk_interval=10): Execute the time-dependent SBM simulation."""
export run_bioreactor_simulation

"""krieger_viscosity(Φ, Φmax): Compute the local fluid viscosity using the Krieger-Dougherty model."""
export krieger_viscosity

"""shear_rate(u): Compute the magnitude of the shear rate tensor for SBM migration."""
export shear_rate

"""particle_flux(u, Φ, gradΦ, params): Compute the particle migration flux combining shear-induced and buoyancy effects."""
export particle_flux

"""navier_stokes_weak_form(u, p, v, q, dΩ, params): Assemble the Navier-Stokes weak form for the bioreactor flow."""
export navier_stokes_weak_form

"""coupled_bioreactor_residual(X, Y, dΩ, params): Compute the monolithic residual for the (u, p, Φ, C) system."""
export coupled_bioreactor_residual

function run_simulation()
    println("Simulation started.")
end

end
