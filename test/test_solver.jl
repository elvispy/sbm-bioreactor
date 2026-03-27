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
    G = TestFESpace(model, reffe_s, conformity=:H1)

    u_fun(x) = VectorValue(0.0, 0.0)
    p_fun(x) = 0.0
    Φ_fun(x) = 0.1
    C_fun(x) = 5.5
    Γ_fun(x) = 0.0

    uh = interpolate_everywhere(u_fun, V)
    ph = interpolate_everywhere(p_fun, Q)
    Φh = interpolate_everywhere(Φ_fun, W)
    Ch = interpolate_everywhere(C_fun, Z)
    Γh = interpolate_everywhere(Γ_fun, G)

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

    x = (uh, ph, Φh, Ch, Γh)
    x_prevs = (uh, ph, Φh, Ch, Γh)
    y = (uh, ph, Φh, Ch, Γh)

    res = coupled_bioreactor_residual(x, (x_prevs,), y, dt, params, 1, t)

    degree = 2
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)

    val = sum(∫(res)dΩ)

    @test typeof(val) == Float64
    @test isfinite(val)

    # Test BDF2 call
    res2 = coupled_bioreactor_residual(x, (x_prevs, x_prevs), y, dt, params, 2, t)
    val2 = sum(∫(res2)dΩ)
    @test typeof(val2) == Float64
    @test isfinite(val2)

    # Any order other than 1 is currently treated as the multi-step branch.
    res3 = coupled_bioreactor_residual(x, (x_prevs, x_prevs), y, dt, params, 99, t)
    val3 = sum(∫(res3)dΩ)
    @test isfinite(val3)

    ablated_params = merge(
        params,
        (
            enable_particle_flux = false,
            freeze_viscosity = true,
            include_convection = false,
            enable_growth_source = false,
            enable_nutrient_reaction = false,
        ),
    )
    res4 = coupled_bioreactor_residual(x, (x_prevs,), y, dt, ablated_params, 1, t)
    val4 = sum(∫(res4)dΩ)
    @test isfinite(val4)

    jac4 = coupled_bioreactor_jacobian(x, (x_prevs,), x, y, dt, merge(ablated_params, (use_explicit_jacobian=true,)), 1, t)
    jac_val4 = sum(∫(jac4)dΩ)
    @test isfinite(jac_val4)

    explicit_params = merge(ablated_params, (use_explicit_jacobian=true,))
    X = MultiFieldFESpace([TrialFESpace(V), TrialFESpace(Q), TrialFESpace(W), TrialFESpace(Z), TrialFESpace(G)])
    Ymf = MultiFieldFESpace([V, Q, W, Z, G])
    op = SBM_Bioreactor.build_bioreactor_operator(X, Ymf, dΩ, (x_prevs,), dt, explicit_params, 1, t)
    @test op isa Gridap.FESpaces.FEOperator

    unsupported_params = merge(params, (use_explicit_jacobian=true,))
    @test_throws ErrorException SBM_Bioreactor.build_bioreactor_operator(X, Ymf, dΩ, (x_prevs,), dt, unsupported_params, 1, t)
end
