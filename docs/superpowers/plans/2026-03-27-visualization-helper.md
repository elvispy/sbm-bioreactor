# Visualization Helper Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `Plots` as a real dependency and make `scripts/visualize.jl` a reusable helper for notebook and standalone HARV visualization.

**Architecture:** Keep visualization logic in a helper script for now, reuse it from the notebook with `include`, and avoid fake VTK-reader behavior.

**Tech Stack:** Julia, Plots.jl, Gridap.jl, Jupyter notebook JSON.

---

## Chunk 1: Visualization Helpers

### Task 1: Replace placeholder visualization code

**Files:**
- Modify: `scripts/visualize.jl`
- Modify: `notebooks/harv_2d_tutorial.ipynb`

- [ ] **Step 1: Add reusable plotting and animation helpers**
- [ ] **Step 2: Update the notebook to call those helpers**
- [ ] **Step 3: Validate notebook JSON**

### Task 2: Declare dependencies and add narrow coverage

**Files:**
- Modify: `Project.toml`
- Modify: `test/test_examples.jl`

- [ ] **Step 1: Add `Plots` as a declared dependency**
- [ ] **Step 2: Add a narrow helper-loading test that only runs when `Plots` is available**
- [ ] **Step 3: Run the narrow example test if the environment supports the new dependency**
