# Unsteady MMS Design

## Goal

Add an unsteady manufactured-solution verification case to the test suite without changing or replacing the existing steady MMS coverage.

## Scope

The new coverage stays inside the existing Julia `test/` suite and reuses the current Gridap-based MMS pattern already present in `test/test_mms.jl`. The steady test remains untouched. The unsteady case is added as a separate test file so both cases coexist cleanly.

## Approach Options

### Option 1: Add a second `@testset` to `test/test_mms.jl`

This keeps all MMS coverage together, but it makes the existing file denser and slightly harder to scan.

### Option 2: Add `test/test_mms_unsteady.jl`

This keeps the new time-dependent logic isolated from the existing steady case while still using the same runtime harness. This is the recommended option because it preserves the current steady file and makes future BDF2 or convergence extensions easier.

### Option 3: Build a reusable MMS helper layer first

This would reduce duplication, but it is unnecessary for a single new test and would expand the change surface without immediate value.

## Recommended Design

Create `test/test_mms_unsteady.jl` and include it from `test/runtests.jl`. The test uses the same five-field mixed formulation and the same manufactured-source construction:

1. Define smooth time-dependent exact fields `u(x,t)`, `p(x,t)`, `Φ(x,t)`, `C(x,t)`, and `Γ(x,t)`.
2. Build the current exact state `x_exact(t)` and previous exact state `x_exact(t-dt)`.
3. Define the internal residual with the actual package residual function.
4. Manufacture the forcing by subtracting the residual evaluated at the exact current state.
5. Verify that the unsteady BDF1 residual is active for a nontrivial previous state and that the manufactured residual vanishes at the exact current-time state.

## Testing

Run the unsteady MMS test directly first, then run the full suite only if the narrow test reveals integration issues that require broader confirmation. Prefer residual-level verification over a full nonlinear solve if the solve is too expensive for a regression test.
