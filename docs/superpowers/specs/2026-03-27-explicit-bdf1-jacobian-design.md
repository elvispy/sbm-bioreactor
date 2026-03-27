# Explicit BDF1 Jacobian Design

## Goal

Replace the remaining autodiff-backed Jacobian path for the real coupled `BDF1` monolithic solve with an explicit weak-form Jacobian, while keeping the current residual and overall solver formulation unchanged.

## Why This Is The Next Step

The strongest performance evidence in the codebase points to Jacobian formation as the main structural bottleneck. The easy-case experiments showed that a hand-written explicit Jacobian can reduce fused residual/Jacobian evaluation from impractical to sub-second on tiny cases. By contrast, solver-method changes such as trust region did not rescue the real coupled model while autodiff still owned the full Jacobian path. The highest-signal next move is therefore to remove the remaining autodiff dependence for `BDF1` rather than introducing new nonlinear strategies or architectural changes.

## Scope

This design covers only the full coupled `BDF1` Jacobian:

- Keep `coupled_bioreactor_residual(...)` as the governing residual.
- Implement the explicit weak-form Jacobian for the real coupled `BDF1` model in package code.
- Use the explicit Jacobian in `build_bioreactor_operator(...)` for supported `BDF1` cases.
- Leave `BDF2` on the current fallback path for now.

This design does not include operator splitting, matrix-free methods, nondimensionalization, or a `BDF2` explicit Jacobian.

## Decomposition

The explicit Jacobian should be split into named helpers that mirror the residual structure:

- `jac_gamma_bdf1(...)`
- `jac_ns_bdf1(...)`
- `jac_continuity_bdf1(...)`
- `jac_phi_bdf1(...)`
- `jac_C_bdf1(...)`

The hard local derivatives should be factored into helper utilities:

- shear-rate directional derivative
- viscosity / viscosity-gradient derivatives
- particle-flux directional derivative

The top-level `coupled_bioreactor_jacobian(...)` should sum these blocks rather than containing a single giant undifferentiated expression.

## Migration Strategy

The migration should be incremental and test-driven.

1. Preserve the current easy-case behavior while refactoring the Jacobian into block helpers.
2. Add the next full-physics `BDF1` block behind failing tests.
3. Compare explicit and autodiff Jacobians on tiny cases after each new block.
4. Expand support until the full coupled `BDF1` path is explicit.
5. Switch the production `BDF1` builder to use the explicit path for real physics.

This keeps the risk localized and makes it easier to identify which block introduces algebra or domain mismatches.

## Testing Strategy

The guardrails should be layered:

1. Tiny-case explicit-vs-autodiff Jacobian comparisons for each migration stage.
2. Package-level explicit operator smoke tests for richer `BDF1` configurations.
3. Constitutive derivative checks in the testbed using symbolic or analytic references for local expressions only.

`Symbolics.jl` is appropriate as a test-only aid for local constitutive expressions such as viscosity derivatives or pieces of the particle flux. It should not be introduced into the production solver path.

## Deferred Alternatives

Operator splitting and matrix-free methods remain valid future alternatives, but they are not the right next step yet. Both require larger architectural changes than the currently best-supported fix. The immediate need is to finish the migration away from autodiff for the real coupled `BDF1` Jacobian first.
