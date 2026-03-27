# HARV Tutorial Notebook Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a reusable 2D HARV example API and a tutorial notebook that visualizes a coarse solver run.

**Architecture:** Put HARV case construction and optional history collection in package code, keep plotting and narrative in the notebook, and verify the new example API with a narrow Julia test.

**Tech Stack:** Julia, Gridap.jl, package test suite, Jupyter `.ipynb` notebook format.

---

## Chunk 1: Reusable Example API

### Task 1: Add a 5-field HARV case builder

**Files:**
- Create: `src/examples.jl`
- Modify: `src/SBM_Bioreactor.jl`
- Test: `test/test_examples.jl`

- [ ] **Step 1: Write a failing example-API test**
- [ ] **Step 2: Implement `build_harv_2d_case` with mapped geometry, spaces, parameters, and metadata**
- [ ] **Step 3: Export the new helper from the package**
- [ ] **Step 4: Re-run the narrow test**

### Task 2: Add optional history collection to the solver loop

**Files:**
- Modify: `src/solver.jl`
- Test: `test/test_examples.jl`

- [ ] **Step 1: Extend `run_bioreactor_simulation` with `collect_history` and configurable VTK prefixing**
- [ ] **Step 2: Update the test to verify zero-step history collection**
- [ ] **Step 3: Re-run the narrow test**

## Chunk 2: Tutorial Artifact

### Task 3: Add the notebook walkthrough

**Files:**
- Create: `notebooks/harv_2d_tutorial.ipynb`

- [ ] **Step 1: Write notebook markdown sections for model, geometry, BCs, and limitations**
- [ ] **Step 2: Add code cells for mesh plotting, initial conditions, coarse run, and animation export**
- [ ] **Step 3: Keep notebook logic thin by calling package helpers**

### Task 4: Wire test runner

**Files:**
- Modify: `test/runtests.jl`

- [ ] **Step 1: Include `test/test_examples.jl` in the test runner**
- [ ] **Step 2: Run the narrow example test**
- [ ] **Step 3: If runtime permits, run a broader package smoke test**
