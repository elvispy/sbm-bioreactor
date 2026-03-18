"""
    scripts/run_harv_2d.jl

Simulation script for a 2D bioreactor using the SBM solver.

This script sets up:
1. A disk-mapped 2D Cartesian grid.
2. Taylor-Hood/Lagrangian finite element spaces.
3. Physical parameters (Krieger-Dougherty viscosity, SBM fluxes).
4. Initial concentration distribution ("floating cells").
5. Execution of the monolithic solver via `run_bioreactor_simulation`.
"""

using Gridap
using SBM_Bioreactor
using LineSearches: BackTracking

# 1. Geometry and Mesh Generation
# The bioreactor domain is a circle of radius 0.05m.
# We map a square [-1,1]^2 to the disk using an elliptical mapping to maintain grid structure.
function square_to_disk(p)
    x, y = p[1], p[2]
    # Map from square grid [-0.05, 0.05] internally to disk coordinates
    X = x / 0.05
    Y = y / 0.05
    x_new = X * sqrt(1.0 - Y^2/2.0)
    y_new = Y * sqrt(1.0 - X^2/2.0)
    return VectorValue(x_new * 0.05, y_new * 0.05)
end

domain = (-0.05, 0.05, -0.05, 0.05)
partition = (10, 10) # 100 elements total
model = CartesianDiscreteModel(domain, partition, map=square_to_disk)

# 2. Finite Element Space Definition
# Using Taylor-Hood inspired elements for stability:
# - Vector-valued lagrangian (2nd order) for velocity
# - Scalar-valued lagrangian (1st order) for pressure/transport
reffe_u = ReferenceFE(lagrangian, VectorValue{2, Float64}, 2)
reffe_p = ReferenceFE(lagrangian, Float64, 1)
reffe_s = ReferenceFE(lagrangian, Float64, 1)

# BCs: Wall rotation u = ω * (y, -x)
ω = 7.5 * 2π / 60.0
u_wall(x) = VectorValue(ω * x[2], -ω * x[1])

# Setup spaces
V = TestFESpace(model, reffe_u, conformity=:H1, dirichlet_tags="boundary")
Q = TestFESpace(model, reffe_p, conformity=:H1, constraint=:zeromean)
W = TestFESpace(model, reffe_s, conformity=:H1)
Z = TestFESpace(model, reffe_s, conformity=:H1)

U = TrialFESpace(V, u_wall)
P = TrialFESpace(Q)
Φ_space = TrialFESpace(W)
C_space = TrialFESpace(Z)

# Multi-field spaces for the monolithic solver
Y = MultiFieldFESpace([V, Q, W, Z])
X = MultiFieldFESpace([U, P, Φ_space, C_space])

# 3. Parameter Definition
degree = 4
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

# Physical parameters as defined in the SBM literature
params = (
    μf = 0.5889,
    Φmax = 0.64,
    a = 5.0e-6,
    ρs = 1000.0,
    ρf = 1050.0,
    g = VectorValue(0.0, -9.81),
    Df = 5.4e-10,
    Φavg = 0.05,
    L = 0.01,
    u_wall = u_wall,
    kc = 1.0e-13,
    ke = 4.2e-6,
    d0 = 3.0e5,
    u0 = u_wall,
    p0 = x -> 0.0,
    Φ0 = x -> (x[2] > 0.02 ? 0.1 : 0.0), # Initial cell cloud setup
    C0 = x -> 5.5
)

# 4. Solver Execution
dt = 0.1
total_time = 10.0
nsteps = Int(total_time / dt)

println("Starting Circular Bioreactor Simulation...")
println("Elements: 100, dt: $dt, Total Steps: $nsteps")

# Run monolithic simulation and save output for visualization
run_bioreactor_simulation(X, Y, dΩ, dt, params, nsteps; write_vtk_interval=10)

println("Simulation finished. Results saved to .vtu files.")
