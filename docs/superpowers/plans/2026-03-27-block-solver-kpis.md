# 2026-03-27 Block Solver KPIs

## Goal

Evaluate a separate blocked linear-solver experiment against the current direct monolithic baseline.

## Baseline

- `64x64` warmed one-step direct solve: `4.23 s`
- `128x128` warmed one-step direct solve: `17.13 s`
- `256x256` warmed one-step direct solve: `89.82 s`

## First milestone

Use the same explicit `BDF1` Jacobian on a blocked `(u,p)` / `(Φ,C,Γ)` layout and compare a block-aware linear solve against the direct baseline.

## KPIs

- Linear solver converges or fails.
- Outer iteration count.
- Inner iteration counts, if exposed.
- Wall-clock linear solve time.
- Residual norm reduction.
- Total nonlinear step time, once the linear probe works.
- Scaling slope from `64x64` to `128x128` to `256x256`.

## Success criteria

- `64x64` converges reliably.
- `64x64` is within roughly `0.5x` to `2x` of direct.
- Scaling is better than direct at larger meshes.

## Partial win

- `64x64` converges but is slower.
- Block structure reveals the limiting block clearly.
- Iteration growth is manageable enough to justify better preconditioners.

## Failure

- No convergence at `64x64`.
- Severe instability or brittle tuning.
- No scaling advantage by `128x128`.
