# Solver Alternatives Note

This note records structural alternatives that remain important even though the current implementation step focuses on an explicit assembled Jacobian for the monolithic solver.

## Deferred Alternatives

### Operator Splitting

The strongest architectural alternative to the monolithic explicit-Jacobian path is a segregated or block-iterative solve. This could reduce nonlinear solve complexity and simplify preconditioning, but it changes the mathematical algorithm and introduces coupling-error and convergence-design questions. It remains a serious future option, not a discarded one.

### Matrix-Free

Matrix-free or JFNK-style methods remain a credible long-term direction for large 3D solves, especially if future work emphasizes Krylov methods, accelerator readiness, or avoiding expensive sparse assembly. The current codebase is not yet ready for that leap because it still lacks a mature preconditioning strategy and a well-controlled linearized operator path.

## Current Choice

The current implementation step prioritizes an explicit assembled Jacobian because it is the lowest-risk structural move that preserves the monolithic formulation while unlocking tractable Newton solves, caching opportunities, and later preconditioning work.
