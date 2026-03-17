module SBM_Bioreactor
using Gridap

include("physics.jl")
include("solver.jl")

export run_simulation
export navier_stokes_weak_form
export krieger_viscosity
export shear_rate
export particle_flux
export coupled_bioreactor_residual

function run_simulation()
    println("Simulation started.")
end

end
