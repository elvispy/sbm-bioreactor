# 2026-03-27 Iterative Block Solver

## Goal

Replace the current direct linear-solve baseline with a first block-aware iterative stress test on the explicit `BDF1` path.

## Why now

- Warm explicit assembly is no longer dominant.
- Large-mesh runtime is solve-dominated.
- The current nonlinear stack still defaults to Gridap's direct `BackslashSolver`.
- The explicit `BDF1` path is now stable enough to benchmark solver structure.

## Baseline evidence

- `64x64`: assembly `0.93 s`, solve `4.25 s`, `2.5M` nnz.
- `128x128`: assembly `2.24 s`, solve `17.13 s`, `10.2M` nnz.
- `256x256`: assembly `8.58 s`, solve `89.82 s`, `40.9M` nnz.

## Chosen next step

Implement a `64x64` linear-solver stress test around the warmed explicit `BDF1` Jacobian:

1. Current direct baseline.
2. A monolithic iterative candidate.
3. A block-aware iterative candidate.

## First block structure

- Flow block: `(u, p)`
- Transport block: `(Φ, C, Γ)`

This is the least-controversial split in the current monolithic system.

## Questions to answer

- Does a monolithic Krylov method beat the direct baseline?
- Does a two-block strategy beat monolithic Krylov?
- Is the flow block saddle-point behavior the limiting factor?
- Are we still dominated by factorization-like work, or by iteration count?

## Immediate implementation target

- Prefer `GridapPETSc` / `GridapSolvers` over a custom solver stack.
- Keep the explicit `BDF1` operator unchanged.
- Change only the linear-solver layer first.
