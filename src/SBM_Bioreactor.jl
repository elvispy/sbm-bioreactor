module SBM_Bioreactor
using Gridap

include("physics.jl")

export run_simulation
export navier_stokes_weak_form
export krieger_viscosity
export shear_rate
export particle_flux

function run_simulation()
    println("Simulation started.")
end

end
