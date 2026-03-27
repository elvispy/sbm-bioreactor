using Test
using Gridap
using SBM_Bioreactor

@testset "Method of Manufactured Solutions (MMS) - Unsteady" begin
    domain = (0, 1, 0, 1)
    partition = (1, 1)
    model = CartesianDiscreteModel(domain, partition)

    Φ_base(x) = sin(pi * x[1]) * sin(pi * x[2])
    C_base(x) = cos(pi * x[1]) * cos(pi * x[2])

    u_ex(x, t) = VectorValue(0.0, 0.0)
    p_ex(x, t) = 0.0
    Φ_ex(x, t) = 0.2 * (1.0 + 0.25 * cos(0.7 * t) * Φ_base(x))
    C_ex(x, t) = 5.5 * (1.0 + 0.05 * sin(0.9 * t) * C_base(x))
    Γ_ex(x, t) = 1.0e-5

    params = (
        μf = 0.5889,
        Φmax = 0.64,
        a = 5.0e-6,
        ρs = 1000.0,
        ρf = 1050.0,
        g = VectorValue(0.0, -9.81),
        Df = 5.4e-10,
        Φavg = 0.2,
        L = 0.01,
        u_wall = x -> VectorValue(0.0, 0.0),
        kc = 1.0e-13,
        ke = 4.2e-6,
        d0 = 3.0e5
    )

    reffe_u = ReferenceFE(lagrangian, VectorValue{2, Float64}, 2)
    reffe_p = ReferenceFE(lagrangian, Float64, 1)
    reffe_s = ReferenceFE(lagrangian, Float64, 1)

    V = TestFESpace(model, reffe_u, conformity = :H1, dirichlet_tags = "boundary")
    Q = TestFESpace(model, reffe_p, conformity = :H1, constraint = :zeromean)
    W = TestFESpace(model, reffe_s, conformity = :H1, dirichlet_tags = "boundary")
    Z = TestFESpace(model, reffe_s, conformity = :H1, dirichlet_tags = "boundary")
    G = TestFESpace(model, reffe_s, conformity = :H1, dirichlet_tags = "boundary")

    t = 0.35
    dt = 0.1
    t_prev = t - dt

    U = TrialFESpace(V, x -> u_ex(x, t))
    P = TrialFESpace(Q)
    Φ_space = TrialFESpace(W, x -> Φ_ex(x, t))
    C_space = TrialFESpace(Z, x -> C_ex(x, t))
    Γ_space = TrialFESpace(G, x -> Γ_ex(x, t))

    Y = MultiFieldFESpace([V, Q, W, Z, G])
    X = MultiFieldFESpace([U, P, Φ_space, C_space, Γ_space])

    Ω = Triangulation(model)
    dΩ = Measure(Ω, 2)

    x_exact = interpolate_everywhere(
        [x -> u_ex(x, t), x -> p_ex(x, t), x -> Φ_ex(x, t), x -> C_ex(x, t), x -> Γ_ex(x, t)],
        X,
    )
    x_prev = interpolate_everywhere(
        [x -> u_ex(x, t_prev), x -> p_ex(x, t_prev), x -> Φ_ex(x, t_prev), x -> C_ex(x, t_prev), x -> Γ_ex(x, t_prev)],
        X,
    )

    y_probe = (
        interpolate_everywhere(x -> u_ex(x, t), V),
        interpolate_everywhere(x -> p_ex(x, t), Q),
        interpolate_everywhere(x -> Φ_ex(x, t), W),
        interpolate_everywhere(x -> C_ex(x, t), Z),
        interpolate_everywhere(x -> Γ_ex(x, t), G),
    )

    res_unsteady = coupled_bioreactor_residual(x_exact, (x_prev,),  y_probe, dt, params, 1, t)
    res_frozen   = coupled_bioreactor_residual(x_exact, (x_exact,), y_probe, dt, params, 1, t)

    unsteady_scalar = sum(∫(res_unsteady)dΩ)
    frozen_scalar = sum(∫(res_frozen)dΩ)

    # The chosen manufactured state must actually activate the time-discrete path.
    @test abs(unsteady_scalar - frozen_scalar) > 1e-6

    res_mms = coupled_bioreactor_residual(x_exact, (x_prev,), y_probe, dt, params, 1, t) - res_unsteady
    mms_scalar = sum(∫(res_mms)dΩ)

    # Manufactured forcing must cancel the full unsteady residual at the exact state.
    @test abs(mms_scalar) < 1e-10
end
