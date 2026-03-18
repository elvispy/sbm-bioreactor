# Design Doc: Bioreactor Simulation (Chao & Das 2015) in Gridap.jl

**Date:** 2026-03-17
**Status:** Approved
**Topic:** Implementing a coupled cell-motion and nutrient-transport mixture model in a rotating HARV bioreactor using Finite Elements.

## 1. Overview
The goal is to reproduce the numerical simulation of a High Aspect Ratio Vessel (HARV) bioreactor. The model treats cells as a continuous volume fraction $\Phi$ suspended in a Newtonian medium, using a mixture continuum framework.

## 2. Mathematical Model

### 2.1 Governing Equations
- **Mixture Momentum (Navier-Stokes):**
  $$\rho \frac{\partial \mathbf{u}}{\partial t} + \rho (\mathbf{u}\cdot\nabla)\mathbf{u} = -\nabla p + \nabla\cdot\left[\mu(\Phi)\left(\nabla\mathbf{u} + \nabla\mathbf{u}^{T}\right)\right] + \rho(\Phi) \mathbf{g}$$
- **Mixture Continuity:** $\nabla \cdot \mathbf{u} = 0$ (assuming nearly matched densities $\rho_s \approx \rho_f$).
- **Cell Transport ($\Phi$):**
  $$\frac{\partial \Phi}{\partial t} + \mathbf{u}\cdot\nabla\Phi = -\nabla \cdot \mathbf{J}_s / \rho_s$$
  Where $\mathbf{J}_s$ includes shear-induced migration fluxes $\nabla(\dot{\gamma}\Phi)$ and hindered settling $f_h \mathbf{u}_{st}\Phi$.
- **Nutrient Transport ($C$):**
  $$\frac{\partial C}{\partial t} + \mathbf{u}\cdot\nabla C = D_f \nabla^2 C + r_c$$

### 2.2 Rheology
- **Krieger-Dougherty Viscosity:** $\mu(\Phi) = \mu_f (1 - \Phi/\Phi_{\max})^{-2.5\Phi_{\max}}$
- **Note on Regularization:** A note will be included in the implementation to handle potential singularities as $\Phi \to \Phi_{\max}$, though no cap is applied initially.

## 3. Architecture & Numerical Strategy

### 3.1 Framework
- **Tool:** Julia with `Gridap.jl`.
- **Approach:** Monolithic (fully coupled) system solved via Newton's method.
- **Differentiability:** Leverages Gridap's Automatic Differentiation (AD) for the Jacobian.

### 3.2 Discretization
- **Element Types (Taylor-Hood):**
    - Velocity ($\mathbf{u}$): P2 (Lagrangian)
    - Pressure ($p$): P1 (Lagrangian)
    - Volume Fraction ($\Phi$): P1 (Lagrangian)
    - Nutrient ($C$): P1 (Lagrangian)
- **Time Stepping:** Implicit Backward Euler for stability in non-linear fluxes.
- **Weak Form:** Integration by parts will be used for the shear-induced migration flux to avoid direct 2nd-order derivatives of velocity.

## 4. Geometry & Boundary Conditions
- **Domain:** 2D Circular cross-section (representing HARV).
- **Velocity BC:** Dirichlet $\mathbf{u} = \omega(y, -x)$ on the outer wall (Rotating wall).
- **Concentration BCs:** Zero-flux (natural Neumann) for $\Phi$ and $C$ at walls.
- **Mesh:** Unstructured triangular mesh generated via `Gmsh` or internal Gridap tools.

## 5. Long-term Extensibility
- **3D Scaling:** Gridap's dimension-independent operators allow for future 3D meshes using the same PDE code.
- **Impeller Potential:** Future support for moving impellers will leverage `GridapEmbedded.jl` (Agard-style immersed boundaries).
- **Solver Performance:** Integration with `GridapPETSc.jl` for large-scale 3D iterative solvers.

## 6. Verification Plan
1. **Fluid Check:** Verify steady-state rotating cylinder velocity profile.
2. **Tracer Check:** Verify advection-diffusion of $\Phi$ with constant viscosity.
3. **Coupling Check:** Verify cell "focusing" orbits resulting from the balance of gravity and shear-induced migration.
