# Explicit BDF1 Jacobian Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the remaining autodiff-backed Jacobian path for the full coupled `BDF1` solver with an explicit weak-form Jacobian.

**Architecture:** Keep the existing monolithic residual intact and add explicit `BDF1` Jacobian blocks that mirror the residual structure. Migrate incrementally, using tiny-case explicit-vs-autodiff tests and local constitutive derivative checks as guardrails.

**Tech Stack:** Julia, Gridap, NLsolve via Gridap, package tests, optional Symbolics.jl in tests only

---

## Chunk 1: Guardrails And Refactor Boundary

### Task 1: Add failing tests for the next full-physics BDF1 expansion

**Files:**
- Modify: `test/test_solver.jl`
- Test: `test/test_solver.jl`

- [ ] **Step 1: Write a failing test for a richer `BDF1` explicit Jacobian case**

Add a tiny-case test that requests `use_explicit_jacobian=true` for a `BDF1` configuration beyond the current easy subset and compares assembly or operator creation against expectations.

- [ ] **Step 2: Run the narrow test to verify failure**

Run: `julia --project=. test/test_solver.jl`
Expected: FAIL in the new explicit-Jacobian coverage test.

- [ ] **Step 3: Refactor the existing easy-case Jacobian into named block helpers**

Move the current easy-case terms into block helpers without expanding physics yet.

- [ ] **Step 4: Run the narrow test suite**

Run: `julia --project=. test/test_solver.jl`
Expected: existing easy-case tests still pass; new richer-case test still fails.

- [ ] **Step 5: Commit**

```bash
git add test/test_solver.jl src/solver.jl
git commit -m "Refactor explicit BDF1 Jacobian into blocks"
```

## Chunk 2: Full-Physics BDF1 Block Expansion

### Task 2: Add constitutive derivative helpers with test-only reference checks

**Files:**
- Modify: `src/solver.jl`
- Modify: `test/test_solver.jl`
- Create or Modify: `test/test_symbolic_jacobian_helpers.jl` (if needed)

- [ ] **Step 1: Write failing derivative-helper tests**

Add focused tests for the next constitutive derivative family, using symbolic or analytic references for local expressions where practical.

- [ ] **Step 2: Run the focused tests to verify failure**

Run the narrow relevant test file(s).
Expected: FAIL in the new constitutive derivative tests.

- [ ] **Step 3: Implement the minimal derivative helper**

Add the helper in `src/solver.jl` without changing unrelated residual logic.

- [ ] **Step 4: Run the focused tests again**

Expected: PASS for the new helper tests.

- [ ] **Step 5: Commit**

```bash
git add src/solver.jl test/test_solver.jl test/test_symbolic_jacobian_helpers.jl
git commit -m "Add constitutive derivative helpers for explicit BDF1 Jacobian"
```

### Task 3: Expand the Jacobian block-by-block

**Files:**
- Modify: `src/solver.jl`
- Modify: `test/test_solver.jl`

- [ ] **Step 1: Add the next failing explicit-vs-autodiff Jacobian comparison**

Start with the next most important block family and add a tiny-case comparison test.

- [ ] **Step 2: Run the narrow test**

Run: `julia --project=. test/test_solver.jl`
Expected: FAIL in the new comparison.

- [ ] **Step 3: Implement the minimal Jacobian block needed**

Update the relevant helper and top-level Jacobian sum.

- [ ] **Step 4: Run the narrow test again**

Expected: PASS for that block comparison.

- [ ] **Step 5: Repeat Steps 1-4 for remaining `BDF1` blocks**

Prioritize variable-viscosity and particle-flux terms last.

- [ ] **Step 6: Commit**

```bash
git add src/solver.jl test/test_solver.jl
git commit -m "Extend explicit BDF1 Jacobian to full coupled physics"
```

## Chunk 3: Production Integration And Verification

### Task 4: Switch production BDF1 operator building to the full explicit path

**Files:**
- Modify: `src/solver.jl`
- Modify: `test/test_solver.jl`
- Modify: `scripts/profile_solver_components.jl` (if needed for verification)

- [ ] **Step 1: Write a failing production-path smoke test**

Add a package-level test that builds or runs a tiny real-physics `BDF1` operator with `use_explicit_jacobian=true`.

- [ ] **Step 2: Run the narrow test to verify failure or unsupported-path behavior**

Run: `julia --project=. test/test_solver.jl`
Expected: FAIL or throw on unsupported explicit path.

- [ ] **Step 3: Update `build_bioreactor_operator(...)` to use the full explicit `BDF1` path**

Keep `BDF2` on the fallback path.

- [ ] **Step 4: Run the narrow tests again**

Expected: PASS for the production-path smoke coverage.

- [ ] **Step 5: Run broader verification**

Run the smallest relevant solver test suite and one tiny profiling comparison.

- [ ] **Step 6: Commit**

```bash
git add src/solver.jl test/test_solver.jl scripts/profile_solver_components.jl
git commit -m "Use explicit BDF1 Jacobian in production solver path"
```

Plan complete and saved to `docs/superpowers/plans/2026-03-27-explicit-bdf1-jacobian.md`. Ready to execute.
