# Bioreactor Simulation Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement a coupled Navier-Stokes, cell transport, and nutrient transport model in Julia using Gridap.jl to simulate the Chao & Das (2015) bioreactor.

**Architecture:** A monolithic (fully coupled) Finite Element solver using Taylor-Hood elements (P2/P1) for velocity-pressure and P1 for concentrations. Non-linear coupling (viscosity and shear-induced migration) is handled via Gridap's automatic differentiation and Newton solver.

**Tech Stack:** Julia, Gridap.jl, Gmsh (via GridapGmsh or internal tools).

---

## Chunk 1: Project Setup & Basic Navier-Stokes

### Task 1: Initialize Project & Define Main Module

**Files:**

- Create: `Project.toml`
- Create: `src/SBM_Bioreactor.jl`
- Test: `test/runtests.jl`

- [ ] **Step 1: Create Project.toml with dependencies**

```toml
[deps]
Gridap = "0.17"
GridapGmsh = "0.5"
LinearAlgebra = ""
```

- [ ] **Step 2: Initialize src/SBM_Bioreactor.jl**

```julia
module SBM_Bioreactor
using Gridap
export run_simulation
end
```

- [ ] **Step 3: Write basic test for module loading**

```julia
using Test
using SBM_Bioreactor
@testset "Module Loading" begin
    @test true
end
```

- [ ] **Step 4: Commit**

```bash
git add Project.toml src/SBM_Bioreactor.jl test/runtests.jl
git commit -m "chore: initialize julia project"
```

### Task 2: Implement Steady-State Navier-Stokes (Rotating Cylinder)

**Files:**

- Create: `src/physics.jl`
- Modify: `src/SBM_Bioreactor.jl`
- Test: `test/test_navier_stokes.jl`

- [ ] **Step 1: Define Navier-Stokes Weak Form in src/physics.jl**

```julia
function navier_stokes_weak_form(u, p, v, q, μ, ρ, f)
    res(u, p, v, q) = (ρ * (u ⋅ ∇(u)) ⋅ v) + (μ * ∇(u) ⊙ ∇(v)) - (p * (∇ ⋅ v)) + (q * (∇ ⋅ u)) - (f ⋅ v)
    return res
end
```

- [ ] **Step 2: Write test for Taylor-Hood stability in test/test_navier_stokes.jl**

```julia
using Gridap
using SBM_Bioreactor
# Verify zero velocity in a static box
```

- [ ] **Step 3: Run test and verify fail**
- [ ] **Step 4: Implement minimal solver in src/SBM_Bioreactor.jl**
- [ ] **Step 5: Run test and verify pass**
- [ ] **Step 6: Commit**

---

## Chunk 2: Rheology & Cell Transport

### Task 3: Implement Krieger-Dougherty Viscosity

**Files:**

- Modify: `src/physics.jl`
- Test: `test/test_rheology.jl`

- [ ] **Step 1: Add viscosity function to src/physics.jl**

```julia
function krieger_viscosity(Φ; μf=0.5889, Φmax=0.63)
    return μf * (1 - Φ/Φmax)^(-2.5*Φmax)
end
```

- [ ] **Step 2: Write test for viscosity values**
- [ ] **Step 3: Commit**

### Task 4: Implement Shear-Induced Migration Fluxes

**Files:**

- Modify: `src/physics.jl`
- Test: `test/test_migration.jl`

- [ ] **Step 1: Define migration weak form in src/physics.jl**
- [ ] **Step 2: Write test for flux divergence**
- [ ] **Step 3: Commit**

---

## Chunk 3: Coupling & Time Integration

### Task 5: Monolithic Coupled Solver

**Files:**

- Create: `src/solver.jl`
- Modify: `src/SBM_Bioreactor.jl`

- [ ] **Step 1: Define the fully coupled residual in src/solver.jl**
- [ ] **Step 2: Implement Backward Euler time-stepping loop**
- [ ] **Step 3: Test coupling with a simple 1D case**
- [ ] **Step 4: Commit**

### Task 6: HARV Bioreactor Case Study (2D)

**Files:**

- Create: `scripts/run_harv_2d.jl`

- [ ] **Step 1: Define HARV geometry and rotating BCs**
- [ ] **Step 2: Run simulation and export VTK**
- [ ] **Step 3: Verify "focusing" orbits visually in ParaView**
- [ ] **Step 4: Commit**

---

## Master Equation & Design Reference

*Consolidated from Chao & Das (2015) and implementation design choices. Append to plan before beginning implementation.*

---

### A. Full Coupled PDE System

The simulation solves a two-phase drift-flux mixture model. The unknowns are: mixture velocity $\mathbf{u}$, pressure $p$, cell volume fraction $\Phi$, and nutrient concentration $C$.

#### A.1 Momentum Equation

$$
\rho \frac{\partial \mathbf{u}}{\partial t}
+ \rho (\mathbf{u}\cdot\nabla)\mathbf{u}
= -\nabla p
- \nabla\cdot\left[\rho c_s (1-c_s)\mathbf{u}_{slip}\mathbf{u}_{slip}\right]
+ \nabla\cdot\left[\mu\left(\nabla\mathbf{u} + \nabla\mathbf{u}^{T}\right)\right]
+ \rho \mathbf{g}
\tag{Eq. 1}
$$

where the quadratic slip-stress term $\rho c_s(1-c_s)\mathbf{u}_{slip}\mathbf{u}_{slip}$ represents momentum exchange between phases. In the dilute limit ($\Phi \to 0$) this term vanishes and Eq. 1 reduces to the Navier–Stokes equation.

The slip velocity is recovered from the particle flux:

$$
\mathbf{u}_{slip} = \frac{\mathbf{J}_s}{\Phi\, \rho_s^\circ\, (1-c_s)}
\tag{from Eq. 9}
$$

#### A.2 Mixture Density and Krieger–Dougherty Viscosity

$$
\rho = (1-\Phi)\rho_f^\circ + \Phi\, \rho_s^\circ
\tag{Eq. 2}
$$

$$
\mu = \mu_f \left(1 - \frac{\Phi}{\Phi_{\max}}\right)^{-2.5\,\Phi_{\max}}
\tag{Eq. 3}
$$

$\Phi_{\max} \approx 0.64$ (random close packing). As $\Phi \to \Phi_{\max}$, $\mu \to \infty$, which acts as a natural kinematic barrier preventing over-packing. This divergence is the basis of **Design Choice 2** below.

#### A.3 Weakly Compressible Continuity Equation *(nonstandard)*

Rather than enforcing $\nabla\cdot\mathbf{u} = 0$, the mixture continuity is:

$$
\nabla\cdot\mathbf{u}
= \frac{\rho_s^\circ - \rho_f^\circ}{\rho_s^\circ\,\rho_f^\circ}\,\nabla\cdot\mathbf{J}_s
\tag{Eq. 8}
$$

When $\rho_s^\circ = \rho_f^\circ$ (density-matched suspension), the right-hand side vanishes and this collapses to the incompressible condition. The right-hand side acts as a volumetric source/sink driven by cell flux divergence.

The original compressible mixture continuity (Eq. 4) that leads to Eq. 8 is:

$$
(\rho_s^\circ-\rho_f^\circ)\left[\nabla\cdot\left(\Phi(1-c_s)\mathbf{u}_{slip}\right)\right]
= \rho_f^\circ\,\nabla\cdot\mathbf{u}
\tag{Eq. 4}
$$

#### A.4 Cell Volume-Fraction Transport

Advection–diffusion form (derived form, Eq. 10):

$$
\frac{\partial \Phi}{\partial t} + \mathbf{u}\cdot\nabla\Phi
= \frac{\rho_s^\circ - \rho_f^\circ}{\rho_s^\circ\,\rho_f^\circ}\,\nabla\cdot\mathbf{J}_s
\tag{Eq. 10}
$$

Note: the right-hand side is identical to the divergence term in Eq. 8, so both equations share the same $\nabla\cdot\mathbf{J}_s$ residual — a useful coupling point in the implementation.

#### A.5 Particle Flux Decomposition

The total diffusive particle flux is:

$$
\mathbf{J}_s = \mathbf{J}_{sc} + \mathbf{J}_{s\mu}
\tag{Eq. 11}
$$

**Shear-induced migration flux** (particle interactions):

$$
\mathbf{J}_{sc} = -a^2\Phi^2\,k_{sc}\,\nabla(\dot{\gamma}\Phi)
\tag{Eq. 12}
$$

**Viscosity-gradient flux** (spatial viscosity variation):

$$
\mathbf{J}_{s\mu} = -a^2\Phi^2\,k_\mu\,\dot{\gamma}\,\nabla(\ln\mu)
\tag{Eq. 13}
$$

with calibrated coefficients $k_{sc} \equiv D_\Phi/a^2 = 0.41$ and $k_\mu \equiv D_\mu/a^2 = 0.62$.

**Full governing flux with hindered settling** (Eq. 14, normalised by $\rho_s^\circ$):

$$
\frac{\mathbf{J}_s}{\rho_s^\circ}
= -\left[0.41\,a^2\,\Phi\,\nabla(\dot{\gamma}\Phi)
       + 0.62\,a^2\,\Phi^2\,\dot{\gamma}\,\nabla(\ln\mu)\right]
- f_h\,\mathbf{u}_{st}\,\Phi
\tag{Eq. 14}
$$

#### A.6 Hindered Settling Velocity *(nonstandard)*

Stokes settling velocity:

$$
\mathbf{u}_{st} = \frac{2a^2(\rho_s^\circ - \rho_f^\circ)}{9\mu}\,\mathbf{g}
\tag{Eq. 17}
$$

Hindered settling correction function:

$$
f_h = \frac{\mu_f\,(1-\Phi_{\text{avg}})}{\mu}
\tag{Eq. 18}
$$

where $\Phi_{\text{avg}}$ is the domain-averaged volume fraction (computed globally at each timestep). This factor $\in (0,1]$ reduces the effective settling velocity as the suspension becomes concentrated — and crucially, as $\Phi \to \Phi_{\max}$, $\mu \to \infty$ so $f_h \to 0$, preventing runaway settling. See **Design Choice 2**.

#### A.7 Shear Rate Scalar

$$
\dot{\gamma}
= \left[\tfrac{1}{2}\,\dot{\boldsymbol{\gamma}}:\dot{\boldsymbol{\gamma}}\right]^{1/2},
\qquad
\dot{\boldsymbol{\gamma}} = \nabla\mathbf{u} + (\nabla\mathbf{u})^T
\tag{Eqs. 19–20}
$$

In 2D Cartesian components:

$$
\dot{\gamma}
= \left[\frac{1}{2}\left(4u_x^2 + 2(u_y+v_x)^2 + 4v_y^2\right)\right]^{1/2}
$$

#### A.8 Nutrient Transport

$$
\frac{\partial C}{\partial t} + \mathbf{u}\cdot\nabla C
= D_f\,\nabla^2 C + r_c
\tag{Eq. 24}
$$

$$
r_c = -\mu_c\cdot d
\tag{Eq. 25}
$$

where $d$ is cell number density (cells m$^{-3}$) and $\mu_c$ is the per-cell consumption rate.

#### A.9 Cell Growth Kinetics

$$
\frac{\partial d}{\partial t}
= k_c\,C\,d_0\,e^{k_e t}
\tag{Eq. 26}
$$

Relationship between cell number density and volume fraction:

$$
\Phi = \frac{\pi}{6}\,d\,a^3
\tag{Eq. 28}
$$

#### A.10 Boundary Conditions

| Boundary | Velocity | Cell flux | Nutrient |
|----------|----------|-----------|---------|
| Rotating wall | $\mathbf{u} = \mathbf{u}_r = \omega(y,-x)$ | $\mathbf{J}_s\cdot\hat{n} = 0$ | $\nabla C\cdot\hat{n} = 0$ |
| All walls | no-slip | no-penetration | no-flux |

Rotational wall velocity: $\mathbf{u}_r = \frac{2\pi \cdot \text{rpm}}{60}(y,-x)$

---

### B. Design Choices & Nonstandard Implementation Notes

#### B.1 Weak Compressibility (Continuity Equation)

**Choice:** Use Eq. 8 ($\nabla\cdot\mathbf{u} = \frac{\rho_s^\circ-\rho_f^\circ}{\rho_s^\circ\rho_f^\circ}\nabla\cdot\mathbf{J}_s$) rather than the standard incompressible $\nabla\cdot\mathbf{u}=0$.

**Motivation:** The two phases have different densities ($\rho_s^\circ = 1000$ kg $m^{-3}$, $\rho_f^\circ = 1050$ kg $m^{-3}$). Local accumulation of cells displaces fluid, producing a net volumetric source term proportional to the particle flux divergence. Ignoring this introduces a mass-conservation error that compounds over long simulation times. In the Gridap implementation this means the pressure Poisson problem has a non-zero RHS that must be assembled from the current $\mathbf{J}_s$ field.

#### B.2 Hindered Settling as Numerical Regularization

**Choice:** Retain the full $f_h = \mu_f(1-\Phi_{\text{avg}})/\mu$ factor in the settling flux rather than setting $f_h=1$.

**Motivation:** As $\Phi \to \Phi_{\max}$, the Krieger–Dougherty viscosity diverges ($\mu\to\infty$), driving $f_h \to 0$ and therefore $f_h\,\mathbf{u}_{st} \to 0$. This prevents the settling flux from driving $\Phi$ above the physical packing limit and acts as a natural stabilizer without requiring explicit flux limiters or clamping. It is both physically correct *and* numerically beneficial.

**Implementation note:** $\Phi_{\text{avg}}$ should be computed as a global integral at the start of each timestep and treated as a frozen scalar for that timestep's assembly.

#### B.3 Iterative Decoupling — "Frozen Flow" Toggle

**Choice:** Structure the solver so that one can toggle between:

- **Fully coupled mode:** $\mathbf{u}$, $\Phi$, $C$ updated simultaneously each step.
- **Frozen-flow mode:** $\mathbf{u}$ held fixed for $N$ sub-steps while only $\Phi$ and $C$ are evolved.

**Motivation:** Cell redistribution occurs on a timescale of minutes to hours while the flow field equilibrates in seconds. During the early "sloshing" phase the velocity field changes rapidly, but once a quasi-steady rotation is established the "slow growth" assumption holds: $\partial\mathbf{u}/\partial t \approx 0$ over many $\Phi$-transport steps. Frozen-flow mode dramatically reduces cost during this phase. The toggle also enables easy sensitivity studies (does the flow field matter at all for concentration gradients?).

**Implementation note:** Expose a `freeze_flow::Bool` parameter in the timestepper configuration. When `true`, skip the momentum-pressure solve and reuse the previous velocity field.

#### B.4 Mesh Refinement Strategy

**Choice:** Use a non-uniform mesh with refinement near the outer rotating wall and near the domain center.

**Motivation:** The key physics occur in two regions:

1. **Wall boundary layer** — velocity gradients are steep; $\dot{\gamma}$ is largest here, which directly drives the shear-induced migration flux $\mathbf{J}_{sc} \propto \nabla(\dot{\gamma}\Phi)$.
2. **Center / axis** — the rotational velocity $\omega r$ approaches zero; this is where cells tend to accumulate and where $\nabla\Phi$ is steepest in the steady state.

Uniform meshes waste degrees of freedom in the bulk. In Gridap, implement via `CartesianDiscreteModel` with a coordinate map, or use a pre-generated mesh with `GmshDiscreteModel`.

#### B.5 Initial Condition — "Floating Cells" Test

**Choice:** Begin with all cells initialised at the top of the domain (buoyant initial condition) rather than a uniform distribution.

**Motivation:** CHO cells in culture medium are *denser* than the medium ($\rho_s^\circ = 1000$ < $\rho_f^\circ = 1050$ kg m$^{-3}$), so they actually *float*. Starting with $\Phi = \Phi_0$ concentrated in the upper half creates the maximum possible gradient in both cell concentration and buoyancy force. This is the most demanding test of the migration flux terms: the rotating flow must overcome buoyancy to redistribute cells, and any errors in the $\mathbf{J}_s$ assembly manifest immediately as non-physical accumulation or diffusion. Use this as the primary verification IC before running biology-relevant cases.

---

### C. Reference Parameters (Table 1, Chao & Das 2015)

| Parameter | Symbol | Value | Unit |
| ----------- | ------ | ----- | ---- |
| Glucose diffusivity | $D_f$ | $5.4\times10^{-10}$ | m² s⁻¹ |
| Initial glucose concentration | $C_0$ | 5.5 | mol m⁻³ |
| Fluid density | $\rho_f^\circ$ | 1050 | kg m⁻³ |
| CHO cell density | $\rho_s^\circ$ | 1000 | kg m⁻³ |
| Cell diameter | $a$ | $5.0\times10^{-6}$ | m |
| Medium viscosity | $\mu_f$ | 0.5889 | Pa s |
| Rotation speed | $\omega$ | 7.5 | rpm |
| Max packing fraction | $\Phi_{\max}$ | 0.64 | — |

**Note on buoyancy sign:** $\rho_s^\circ < \rho_f^\circ$, so the net buoyancy force on cells is *upward*. The gravitational term in Eq. 14 therefore drives cells *toward the top* of the domain in the absence of flow. The rotating flow creates the shear-induced migration that counteracts this and distributes cells toward the center — the central result of the Chao & Das paper.

#### B.6 Depth-Averaged Friction (Hele-Shaw term)

**Choice:** Add a body force term $\mathbf{f}_{drag} = \frac{4\mu_f(\mathbf{u}_r - \mathbf{u})}{L^2}$ to the momentum equation (re-derived from Eq. 22).

**Motivation:** The HARV is a thin disk, not an infinite cylinder. In a 2D simulation, the drag from the front and back "windows" is lost. The authors include Eq. 22 to account for this. Without this term, the fluid would eventually reach a rigid-body rotation; with it, the "depth-averaged" velocity $\mathbf{u}$ is constantly being pulled toward the wall velocity $\mathbf{u}_r$ by the stationary/rotating faces of the disk. This is what generates the specific radial pressure gradients needed for the "focusing" orbits.

**Implementation note:** Add `(4.0 * μf / L^2) * (u_wall - u)` as a source term in the `navier_stokes_weak_form`.

#### B.7 Explicit Time-Dependent Growth Kinetics

**Choice:** Use Eq. 26: $\frac{\partial d}{\partial t} = k_c \cdot C \cdot d_0 \cdot e^{k_e t}$.

**Motivation:** This is a "global clock" growth model rather than a standard "local density" model ($\dot{d} = \mu d$). It assumes that the growth potential of the population is a function of time from inoculation ($t=0$), only modulated locally by nutrient availability $C$. While unusual for modern bioprocessing, it is essential for reproducing the Chao 2015 results.

---

### Revolutionary Improvements & Future Frontiers

Beyond reproducing the 2015 paper, the following "high-signal" improvements would transition this model from a research tool to a revolutionary **Digital Twin** framework for bioprocessing:

1. **Mechanobiological Feedback Loops:**
   - *The Idea:* Couple the local shear rate $\dot{\gamma}$ and mechanical pressure $p$ directly into the growth kinetics $k_e$ and consumption $r_c$.
   - *Revolutionary Impact:* This moves from "cells in a fluid" to "cells *responding* to a fluid." It allows for predicting cell death or differentiation in high-shear zones, which is the primary bottleneck in scaling up industrial bioreactors for cell therapies, cultivated meat, and high-value biologics.

2. **Differentiable Simulation & Adjoint-based Design:**
   - *The Idea:* Ensure the entire Gridap solver chain is compatible with Julia's Automatic Differentiation (AD) ecosystem (e.g., ForwardDiff.jl).
   - *Revolutionary Impact:* This enables **Automated Process Optimization**. Instead of manual trial-and-error, we can use the gradient of "Total Yield" to mathematically derive the optimal rotation schedule or "perfect" impeller shape for a specific cell type's shear tolerance.

3. **Resolved Gas-Liquid Interface (Oxygenation):**
   - *The Idea:* Integrate a Level-Set or VoF solver to track the free surface or discrete air bubbles.
   - *Revolutionary Impact:* Oxygen transfer is the universal limiting factor in large-scale bioprocessing. A model that captures the physical "gulping" of air at the surface or bubble-cell collisions would be a gold-standard tool for bioreactor design across the entire pharma and biotech industry.

4. **Hybrid Continuum-Agent Coupling (HCA):**
   - *The Idea:* Use Agent-Based Modeling (ABM) for individual cells in low-density regions and the $\Phi$ mixture model for high-density regions.
   - *Revolutionary Impact:* This captures the stochasticity of early-stage growth (where individual "founder" cells matter) while maintaining the computational efficiency of the continuum model for mass production.
