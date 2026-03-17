using Test
using SBM_Bioreactor

@testset "Rheology (Krieger-Dougherty)" begin
    μf = 0.5889
    Φmax = 0.64
    
    # Test zero volume fraction recovers fluid viscosity
    @test krieger_viscosity(0.0; μf=μf, Φmax=Φmax) ≈ μf
    
    # Test viscosity increases with volume fraction
    μ1 = krieger_viscosity(0.1; μf=μf, Φmax=Φmax)
    μ2 = krieger_viscosity(0.3; μf=μf, Φmax=Φmax)
    @test μ1 > μf
    @test μ2 > μ1
    
    # Test it approaches infinity near maximum packing
    μ_near_max = krieger_viscosity(0.639; μf=μf, Φmax=Φmax)
    @test μ_near_max > 100 * μf
end
