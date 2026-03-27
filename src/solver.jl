using Gridap
using LineSearches: BackTracking

"""
    coupled_bioreactor_residual(x, x_prevs, y, dt, params, order=1, t=0.0)

Compute the residual of the monolithic 5-variable system for the SBM bioreactor.

# 5-Variable Formulation (u, p, Φ, C, Γ)
- `u`: Velocity field (Mixture momentum).
- `p`: Pressure field (Continuity/Incompressibility).
- `Φ`: Particle volume fraction (Particle transport).
- `C`: Nutrient concentration (Chemical transport).
- `Γ`: Shear rate (Projected auxiliary variable for smoothing).

# Physical Meaning of Residuals
1. `res_ns`: Navier-Stokes momentum balance with Hele-Shaw drag.
2. `res_continuity_rhs`: Modified continuity equation. In SBM, divergence of velocity 
   is non-zero if there is a net particle flux due to density differences (ρs != ρf).
3. `res_phi`: Particle conservation law including shear-induced migration and 
   cell growth (source term).
4. `res_C`: Nutrient conservation with advection, diffusion, and consumption.
5. `res_gamma`: L2-projection of the shear rate. Treating Γ as a primary variable 
   smoothes the shear-gradient terms in the particle flux.
"""
function coupled_bioreactor_residual(x, x_prevs, y, dt, params, order=1, t=0.0)
    # Unpack trial and test functions
    u, p, Φ, C, Γ = x
    v, q, w, z, v_γ = y
    
    # Time discretization using BDF1 (Backward Euler) or BDF2
    if order == 1
        u_n, p_n, Φ_n, C_n, Γ_n = x_prevs[1]
        u_dot = (u - u_n) / dt
        Φ_dot = (Φ - Φ_n) / dt
        C_dot = (C - C_n) / dt
    else
        u_n, p_n, Φ_n, C_n, Γ_n = x_prevs[1]
        u_nn, p_nn, Φ_nn, C_nn, Γ_nn = x_prevs[2]
        u_dot = (3.0*u - 4.0*u_n + u_nn) / (2.0*dt)
        Φ_dot = (3.0*Φ - 4.0*Φ_n + Φ_nn) / (2.0*dt)
        C_dot = (3.0*C - 4.0*C_n + C_nn) / (2.0*dt)
    end

    # Parameters from the configuration
    μf = params.μf
    Φmax = params.Φmax
    a = params.a
    ρs = params.ρs
    ρf = params.ρf
    g = params.g
    Df = params.Df
    Φavg = params.Φavg
    L = params.L
    u_wall = params.u_wall
    kc = params.kc
    ke = params.ke
    d0 = params.d0
    freeze_viscosity = hasproperty(params, :freeze_viscosity) ? params.freeze_viscosity : false
    enable_particle_flux = hasproperty(params, :enable_particle_flux) ? params.enable_particle_flux : true
    enable_growth_source = hasproperty(params, :enable_growth_source) ? params.enable_growth_source : true
    enable_nutrient_reaction = hasproperty(params, :enable_nutrient_reaction) ? params.enable_nutrient_reaction : true
    include_convection = hasproperty(params, :include_convection) ? params.include_convection : true
    
    # Constitutive Relations: Mixture density and Krieger-Dougherty viscosity
    ρ = (1.0 - Φ) * ρf + Φ * ρs
    if freeze_viscosity
        μ = μf
        ∇μ = 0.0 * ∇(Φ)
    else
        visc_op(phi) = krieger_viscosity(phi; μf=μf, Φmax=Φmax)
        μ = visc_op ∘ Φ
        
        # Analytical gradient of viscosity with respect to Φ (for the migration terms)
        # dμ/dΦ = μf * (-2.5*Φmax) * (1-Φ/Φmax)^(-2.5*Φmax-1) * (-1/Φmax)
        dμ_dΦ_op(phi) = 2.5 * μf * (1.0 - phi/Φmax)^(-2.5*Φmax - 1.0)
        ∇μ = (dμ_dΦ_op ∘ Φ) * ∇(Φ)
    end
    
    # 1. Shear Rate Projection: Smooth Γ to calculate ∇Γ in the migration flux
    res_gamma = v_γ * (Γ - shear_rate(u))

    # 2. Momentum Balance: Navier Stokes + Hele-Shaw depth-averaged friction
    drag_coeff = 4.0 * μf / (L^2)
    res_ns = (ρ * u_dot ⋅ v) + navier_stokes_weak_form(u, p, v, q, μ, ρ, g, u_wall, drag_coeff; include_convection=include_convection)
    
    # 3. Modified Continuity: ∇⋅u = - ∇⋅(J * (1/ρs - 1/ρf))
    # This accounts for volume changes when particles migrate in a variable density mixture.
    flux = enable_particle_flux ? particle_flux(u, Φ, ∇(Φ), μ, ∇μ, a, ρs, ρf, μf, Φavg, g, Γ, ∇(Γ)) : (0.0 * ∇(Φ))
    res_continuity_rhs = ∇(q) ⋅ (flux * ((ρs - ρf) / (ρs * ρf)))
    
    # 4. Particle Transport: ∂Φ/∂t + u⋅∇Φ = -∇⋅J + Source
    # Source term models cell proliferation (B.7)
    source_phi = enable_growth_source ? (π/6.0 * a^3) * kc * C * d0 * exp(ke * t) : 0.0
    res_phi = (w * Φ_dot) + (w * (u ⋅ ∇(Φ))) + (∇(w) ⋅ (flux * ((ρs - ρf) / (ρs * ρf)))) - (w * source_phi)
    
    # 5. Nutrient Transport: Advection-Diffusion-Reaction
    # Consumption rate rc is proportional to cell concentration (Φ / volume_of_one_cell)
    rc = enable_nutrient_reaction ? -kc * (Φ / (π/6.0 * a^3)) : 0.0
    res_C = (z * C_dot) + (z * (u ⋅ ∇(C))) + (Df * ∇(z) ⊙ ∇(C)) - (z * rc)
    
    return res_ns + res_continuity_rhs + res_phi + res_C + res_gamma
end

"""
    run_bioreactor_simulation(X, Y, dΩ, dt, params, nsteps; write_vtk_interval=1, output_prefix="results", collect_history=false, profile_steps=false)

Execute the time-stepping loop for the bioreactor simulation using a Newton solver.

# Arguments
- `X`, `Y`: Trial and Test MultiFieldFESpaces.
- `dΩ`: Integration measure.
- `dt`: Time step size.
- `params`: NamedTuple of physical and numerical parameters.
- `nsteps`: Number of time steps.
- `write_vtk_interval`: Frequency of VTK output.
"""
function run_bioreactor_simulation(
    X,
    Y,
    dΩ,
    dt,
    params,
    nsteps;
    write_vtk_interval=1,
    output_prefix="results",
    collect_history=false,
    profile_steps=false,
)
    # Initial state interpolation
    x_n = nothing
    initial_setup_time = @elapsed begin
        x_n = interpolate_everywhere([params.u0, params.p0, params.Φ0, params.C0, params.Γ0], X)
    end
    if profile_steps
        println("profile: initial_setup_time=$(round(initial_setup_time, digits=3)) s")
    end
    x_nn = x_n # For BDF2, first step fallback to BDF1 logic
    
    # Newton-Raphson solver with BackTracking line search for stability
    nls = NLSolver(show_trace=true, method=:newton, linesearch=BackTracking())
    solver = FESolver(nls)
    
    xh = x_n
    history = collect_history ? Any[x_n] : nothing
    times = collect_history ? Float64[0.0] : nothing
    step_profiles = profile_steps ? NamedTuple[] : nothing
    
    for step in 1:nsteps
        t = step * dt
        println("Step: $step, Time: $t")
        
        # Use BDF1 for the first step, BDF2 for subsequent steps
        order = step == 1 ? 1 : 2
        x_prevs = order == 1 ? (x_n,) : (x_n, x_nn)
        
        # Define the residual on the triangulation
        op = nothing
        op_build_time = @elapsed begin
            res(x, y) = ∫( coupled_bioreactor_residual(x, x_prevs, y, dt, params, order, t) )dΩ
            op = FEOperator(res, X, Y)
        end
        if profile_steps
            println("profile: step=$(step) operator_build_time=$(round(op_build_time, digits=3)) s")
        end
        
        # Solve the nonlinear system
        solve_result = nothing
        solve_time = @elapsed begin
            solve_result = solve!(xh, solver, op)
        end
        if profile_steps
            println("profile: step=$(step) solve_time=$(round(solve_time, digits=3)) s")
        end
        xh, _ = solve_result
        
        # Update time-history
        x_nn = x_n
        x_n = xh
        if collect_history
            push!(history, xh)
            push!(times, t)
        end
        
        # Diagnostic output
        vtk_time = 0.0
        if write_vtk_interval > 0 && step % write_vtk_interval == 0
            vtk_time = @elapsed begin
                writevtk(get_triangulation(dΩ), "$(output_prefix)_$step",
                         cellfields=["u"=>xh[1], "p"=>xh[2], "phi"=>xh[3], "C"=>xh[4], "gamma"=>xh[5]])
            end
        end
        if profile_steps && vtk_time > 0
            println("profile: step=$(step) vtk_time=$(round(vtk_time, digits=3)) s")
        end

        if profile_steps
            push!(step_profiles, (
                step = step,
                time = t,
                order = order,
                operator_build_time = op_build_time,
                solve_time = solve_time,
                vtk_time = vtk_time,
            ))
        end
    end
    
    if collect_history
        if profile_steps
            return (
                final_state = xh,
                history = history,
                times = times,
                profile = (
                    initial_setup_time = initial_setup_time,
                    steps = step_profiles,
                ),
            )
        end
        return (final_state=xh, history=history, times=times)
    end
    if profile_steps
        return (
            final_state = xh,
            profile = (
                initial_setup_time = initial_setup_time,
                steps = step_profiles,
            ),
        )
    end
    return xh
end
