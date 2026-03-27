# Visualization Helper Design

## Goal

Make `Plots.jl` a declared project dependency and turn `scripts/visualize.jl` into a useful visualization helper for both the HARV tutorial notebook and standalone script usage.

## Scope

This change adds real plotting and animation helpers for in-memory HARV tutorial workflows. It does not implement generic VTK postprocessing in this pass, because the current script only contains placeholders and no verified reader path.

## Approach Options

### Option 1: Recommended, reusable helper script plus notebook reuse

Keep visualization functions in `scripts/visualize.jl` and call them from the notebook via `include`. This provides immediate utility with minimal restructuring.

### Option 2: Move visualization fully into `src/`

This would be cleaner long-term, but it promotes plotting into the public package surface earlier than needed.

### Option 3: Leave notebook plotting inline and only patch the script

This would preserve duplication and guarantee future drift.

## Recommended Design

1. Add `Plots` to `Project.toml`.
2. Replace placeholder VTK comments in `scripts/visualize.jl` with real helpers for:
   - mapped HARV mesh visualization,
   - scalar-field sampling,
   - initial-condition plots,
   - FE scalar snapshots,
   - GIF/MP4 animation from in-memory history,
   - a small standalone demo function based on `build_harv_2d_case`.
3. Update the tutorial notebook to reuse those helpers instead of duplicating plotting code.

## Testing

1. Validate the notebook JSON after edits.
2. Add a narrow test covering helper loading and at least one non-run plotting path.
3. If runtime permits, instantiate the plotting dependency and run the narrow test. If dependency installation is blocked, report that clearly.
