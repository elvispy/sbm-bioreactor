# Unsteady MMS Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a separate unsteady MMS regression test that coexists with the existing steady MMS test.

**Architecture:** Reuse the current Gridap MMS harness and isolate the new time-dependent manufactured solution in its own test file. Wire the new file into the normal test runner without changing solver or physics code.

**Tech Stack:** Julia, Gridap.jl, Test stdlib, existing `SBM_Bioreactor` residual and nonlinear solve path.

---

## Chunk 1: Test Surface

### Task 1: Add the unsteady MMS test file

**Files:**
- Create: `test/test_mms_unsteady.jl`
- Modify: `test/runtests.jl`
- Reference: `test/test_mms.jl`

- [ ] **Step 1: Write the failing unsteady MMS test**

Add a new testset with time-dependent manufactured fields, previous-state construction, residual subtraction, and residual-level assertions proving the unsteady BDF1 path is exercised.

- [ ] **Step 2: Run the narrow unsteady MMS test**

Run: `julia --project=. test/test_mms_unsteady.jl`
Expected: fail if the harness or manufactured-state construction is incomplete or incorrect.

- [ ] **Step 3: Adjust the test to match the existing Gridap/solver conventions**

Keep the steady MMS untouched. Use the same FE spaces and residual shape already used in `test/test_mms.jl`.

- [ ] **Step 4: Wire the new test into the main runner**

Include `test/test_mms_unsteady.jl` from `test/runtests.jl`.

- [ ] **Step 5: Re-run the narrow MMS verification**

Run: `julia --project=. test/runtests.jl`
Expected: the new unsteady MMS test and existing steady MMS test both pass.
