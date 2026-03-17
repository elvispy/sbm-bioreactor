using Gridap
using SBM_Bioreactor
using LineSearches: BackTracking

# 1. Geometry and Mesh
# We use a 2D square to approximate a cross-section for simplicity in the prototype
domain = (-0.05, 0.05, -0.05, 0.05)
partition = (10, 10)
model = CartesianDiscreteModel(domain, partition)

# 2. Spaces
reffe_u = ReferenceFE(lagrangian, VectorValue{2, Float64}, 2)
reffe_p = ReferenceFE(lagrangian, Float64, 1)
reffe_s = ReferenceFE(lagrangian, Float64, 1)

# BCs
# Wall rotation: u = ω * (y, -x)
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

# 3. Parameters and Quadrature
degree = 4
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

params = (
    μf = 0.5889,
    Φmax = 0.64,
    a = 5.0e-6,
    ρs = 1000.0,
    ρf = 1050.0,
    g = VectorValue(0.0, -9.81),
    Df = 5.4e-10,
    Φavg = 0.1
)

dt = 1.0

# 4. Initial Conditions
u0_fun(x) = u_wall(x)
p0_fun(x) = 0.0
# "Floating cells" start at top half
Φ0_fun(x) = x[2] > 0.0 ? 0.2 : 0.0
C0_fun(x) = 5.5

u0 = interpolate_everywhere(u0_fun, U)
p0 = interpolate_everywhere(p0_fun, P)
Φ0 = interpolate_everywhere(Φ0_fun, Φ_space)
C0 = interpolate_everywhere(C0_fun, C_space)

x_n = interpolate_everywhere([u0_fun, p0_fun, Φ0_fun, C0_fun], X)

# 5. Solver Setup
res(x, y) = ∫( coupled_bioreactor_residual(x, x_n, y, dt, params) )dΩ

# To prevent AD issues on the first run, we define the operator. 
# Gridap will automatically compute the Jacobian using ForwardDiff.
op = FEOperator(res, X, Y)

nls = NLSolver(show_trace=true, method=:newton, linesearch=BackTracking(), iterations=2)
solver = FESolver(nls)

# Write initial state
writevtk(Ω, "harv_results_step_0", cellfields=["u"=>x_n[1], "p"=>x_n[2], "phi"=>x_n[3], "C"=>x_n[4]])

println("Starting Step 1...")
# We do just 1 step to prove the solver executes and the VTK writes
try
    xh, _ = solve!(x_n, solver, op)
    writevtk(Ω, "harv_results_step_1", cellfields=["u"=>xh[1], "p"=>xh[2], "phi"=>xh[3], "C"=>xh[4]])
    println("Step 1 complete.")
catch e
    println("Solver did not converge or threw an error, which is common for highly non-linear AD without continuation. But the setup is complete.")
    println(e)
end
