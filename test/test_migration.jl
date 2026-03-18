using Test
using Gridap
using SBM_Bioreactor

@testset "Shear-Induced Migration Fluxes" begin
    # Test shear rate function
    # Let's create a simple velocity field u(x,y) = (y, 0)
    # Then ∇u = [0 1; 0 0]
    # ε(u) = 1/2 (∇u + ∇u^T) = [0 0.5; 0.5 0]
    # ε ⊙ ε = 0.25 + 0.25 = 0.5
    # γ̇ = sqrt(2 * 0.5) = 1.0

    domain = (0, 1, 0, 1)
    partition = (4, 4)
    model = CartesianDiscreteModel(domain, partition)
    
    reffe_u = ReferenceFE(lagrangian, VectorValue{2, Float64}, 2)
    V = TestFESpace(model, reffe_u, conformity=:H1)
    
    # Velocity field u(x) = (x[2], 0.0) -> simple shear
    u_fun(x) = VectorValue(x[2], 0.0)
    uh = interpolate_everywhere(u_fun, V)
    
    degree = 2
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)
    
    # Check that shear_rate integrates to 1.0 roughly
    # Area is 1.0, so ∫ γ̇ dΩ should be 1.0
    val_gamma = sum(∫(shear_rate(uh))dΩ)
    
    @test isapprox(val_gamma, 1.0, atol=1e-4)

    # Test particle flux components (we just test they run without throwing for now)
    # We need a Φ field and μ field
    reffe_s = ReferenceFE(lagrangian, Float64, 1)
    Q = TestFESpace(model, reffe_s, conformity=:H1)
    
    Φ_fun(x) = 0.1
    Φh = interpolate_everywhere(Φ_fun, Q)
    
    μ_fun(x) = 1.0
    μh = interpolate_everywhere(μ_fun, Q)

    a = 5.0e-6
    ρs = 1000.0
    ρf = 1050.0
    μf = 0.5889
    Φavg = 0.1
    g = VectorValue(0.0, -9.81)
    
    Γ_fun(x) = 1.0
    Γh = interpolate_everywhere(Γ_fun, Q)

    flux = particle_flux(uh, Φh, ∇(Φh), μh, ∇(μh), a, ρs, ρf, μf, Φavg, g, Γh, ∇(Γh))
    
    # Integrate the flux over the domain to see it evaluates
    flux_int = sum(∫(flux)dΩ)
    
    @test typeof(flux_int) <: VectorValue{2, Float64}
end
