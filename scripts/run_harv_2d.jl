using Gridap
using SBM_Bioreactor
using LineSearches: BackTracking

# 1. Geometry and Mesh
domain = (-0.05, 0.05, -0.05, 0.05)
partition = (4, 4)
model = CartesianDiscreteModel(domain, partition)

# 2. Spaces
reffe_u = ReferenceFE(lagrangian, VectorValue{2, Float64}, 2)
reffe_p = ReferenceFE(lagrangian, Float64, 1)
reffe_s = ReferenceFE(lagrangian, Float64, 1)

# BCs
ω = 7.5 * 2π / 60.0
u_wall(x) = VectorValue(ω * x[2], -ω * x[1])

V = TestFESpace(model, reffe_u, conformity=:H1, dirichlet_tags="boundary")
Q = TestFESpace(model, reffe_p, conformity=:H1, constraint=:zeromean)
W = TestFESpace(model, reffe_s, conformity=:H1)
Z = TestFESpace(model, reffe_s, conformity=:H1)

U = TrialFESpace(V, u_wall)
P = TrialFESpace(Q)
Φ_space = TrialFESpace(W)
C_space = TrialFESpace(Z)

Y = MultiFieldFESpace([V, Q, W, Z])
X = MultiFieldFESpace([U, P, Φ_space, C_space])

# 3. Parameters
degree = 4
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

# Growth parameters
params = (
    μf = 0.5889,
    Φmax = 0.64,
    a = 5.0e-6,
    ρs = 1000.0,
    ρf = 1050.0,
    g = VectorValue(0.0, -9.81),
    Df = 5.4e-10,
    Φavg = 0.1,
    u0 = u_wall,
    p0 = x -> 0.0,
    Φ0 = x -> (x[2] > 0.0 ? 0.2 : 0.0),
    C0 = x -> 5.5
)

dt = 0.1
nsteps = 2

println("Starting Simulation with BDF2...")
run_bioreactor_simulation(X, Y, dΩ, dt, params, nsteps; write_vtk_interval=1)
println("Simulation finished.")
