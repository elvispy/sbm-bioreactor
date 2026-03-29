using Test
using Gridap
using SBM_Bioreactor

@testset "Coupled Solver Fast Checks" begin
    domain = (0, 1, 0, 1)
    partition = (2, 2)
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
    Γ_fun(x) = 1.0e-5

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
    x_prev = x
    y = x

    Ω = Triangulation(model)
    dΩ = Measure(Ω, 2)

    res = coupled_bioreactor_residual(x, (x_prev,), y, dt, params, 1, t)
    @test isfinite(sum(∫(res)dΩ))

    ablated_params = merge(
        params,
        (
            enable_particle_flux = false,
            freeze_viscosity = true,
            include_convection = false,
            enable_growth_source = false,
            enable_nutrient_reaction = false,
            use_explicit_jacobian = true,
        ),
    )
    jac = coupled_bioreactor_jacobian(x, (x_prev,), x, y, dt, ablated_params, 1, t)
    @test isfinite(sum(∫(jac)dΩ))

    case = build_harv_2d_case(partition=(1, 1), dt=0.25, total_time=0.25, degree=2)
    x0 = interpolate_everywhere([case.params.u0, case.params.p0, case.params.Φ0, case.params.C0, case.params.Γ0], case.X)
    x0_prevs = (x0,)

    full_explicit_params = merge(case.params, (use_explicit_jacobian = true,))
    op = SBM_Bioreactor.build_bioreactor_operator(case.X, case.Y, case.dΩ, x0_prevs, case.metadata.dt, full_explicit_params, 1, case.metadata.dt)
    @test op isa Gridap.FESpaces.FEOperator

    @test SBM_Bioreactor._build_block_linear_solver(transport_kind=:lu, outer_kind=:gmres) isa Gridap.Algebra.LinearSolver
    @test SBM_Bioreactor._build_block_linear_solver(transport_kind=:gmres, outer_kind=:gmres) isa Gridap.Algebra.LinearSolver
    @test_throws ErrorException SBM_Bioreactor._build_block_linear_solver(transport_kind=:does_not_exist, outer_kind=:gmres)
end
