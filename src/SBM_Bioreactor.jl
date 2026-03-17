module SBM_Bioreactor
using Gridap

include("physics.jl")

export run_simulation
export navier_stokes_weak_form
export krieger_viscosity

function run_simulation()
    println("Simulation started.")
end

end
