using Test
using Gridap
using SBM_Bioreactor

@testset "Coupled Solver Residual" begin
    domain = (0, 1, 0, 1)
    partition = (4, 4)
    model = CartesianDiscreteModel(domain, partition)
    
    reffe_u = ReferenceFE(lagrangian, VectorValue{2, Float64}, 2)
    reffe_p = ReferenceFE(lagrangian, Float64, 1)
    reffe_s = ReferenceFE(lagrangian, Float64, 1)
    
    V = TestFESpace(model, reffe_u, conformity=:H1)
    Q = TestFESpace(model, reffe_p, conformity=:H1)
    W = TestFESpace(model, reffe_s, conformity=:H1)
    Z = TestFESpace(model, reffe_s, conformity=:H1)
    
    u_fun(x) = VectorValue(0.0, 0.0)
    p_fun(x) = 0.0
    Φ_fun(x) = 0.1
    C_fun(x) = 5.5
    
    uh = interpolate_everywhere(u_fun, V)
    ph = interpolate_everywhere(p_fun, Q)
    Φh = interpolate_everywhere(Φ_fun, W)
    Ch = interpolate_everywhere(C_fun, Z)
    
    params = (
        μf = 0.5889,
        Φmax = 0.64,
        a = 5.0e-6,
        ρs = 1000.0,
        ρf = 1050.0,
        g = VectorValue(0.0, -9.81),
        Df = 5.4e-10,
        Φavg = 0.1,
        L = 0.01,
        u_wall = x -> VectorValue(0.0, 0.0),
        kc = 1.0e-13,
        ke = 4.2e-6,
        d0 = 3.0e5
    )
    
    dt = 1.0
    t = 0.0
    
    x = (uh, ph, Φh, Ch)
    x_prevs = (uh, ph, Φh, Ch)
    y = (uh, ph, Φh, Ch) # using functions from the space as test functions for simple evaluation
    
    res = coupled_bioreactor_residual(x, (x_prevs,), y, dt, params, 1, t)
    
    degree = 2
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)
    
    val = sum(∫(res)dΩ)
    
    @test typeof(val) == Float64
    
    # Test BDF2 call
    res2 = coupled_bioreactor_residual(x, (x_prevs, x_prevs), y, dt, params, 2)
    val2 = sum(∫(res2)dΩ)
    @test typeof(val2) == Float64
end
