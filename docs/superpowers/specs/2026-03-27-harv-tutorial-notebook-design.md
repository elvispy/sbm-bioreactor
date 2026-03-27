# HARV Tutorial Notebook Design

## Goal

Add a first-user tutorial notebook for the 2D HARV example that explains the solver, shows the mesh and initial conditions, runs a coarse simulation, and produces visual outputs including an animation.

## Scope

The notebook is primarily pedagogical, not a formal benchmark or regression harness. It should remain thin by calling reusable package functions for HARV setup and solver execution. The stale `scripts/run_harv_2d.jl` file is left unchanged in this pass.

## Approach Options

### Option 1: Notebook-first with all logic inline

This would be fast to author, but it would duplicate solver setup and drift quickly as the code evolves.

### Option 2: Recommended, thin notebook over reusable example API

Add package helpers for constructing a valid 5-field HARV case and for optionally collecting time history during the solve. Keep plotting and narrative in the notebook.

### Option 3: Full docs pipeline with generated notebooks

This is maintainable at larger scale but premature for the current repository.

## Recommended Design

1. Add reusable HARV example helpers to package code.
2. Extend `run_bioreactor_simulation` with optional history collection and configurable output prefixing.
3. Add a small test that verifies the HARV example builder returns a consistent 5-field configuration and that history collection works for a zero-step run.
4. Add a tutorial notebook under `notebooks/` with:
   - model overview,
   - LaTeX equations,
   - mapped-disk geometry and mesh visualization,
   - initial-condition plots,
   - one coarse run,
   - scalar field snapshots,
   - GIF/MP4 export cell for a coarse animation.

## Design Constraints

- Do not rely on the stale `scripts/run_harv_2d.jl` path.
- Keep plotting/tooling dependencies inside the notebook where possible.
- Make the notebook honest about being a coarse tutorial run rather than a production configuration.

## Testing

1. Run a targeted Julia test file for the example API.
2. Verify package loading still works after the new exports and helper include.
3. Notebook execution is not required for this implementation pass; instead, provide a notebook that is internally coherent and points out any needed optional packages.
