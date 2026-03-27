# Diary

## 2026-03-27

- Autodiff Jacobian isolated as dominant bottleneck.
- Easy explicit Jacobian proved huge assembly speedup.
- Remote crashes traced to Julia version mismatch.
- `mac-zerotier` downgraded from `1.12.5` to `1.12.1`.
- Full real-physics `BDF1` explicit operator now builds.
- Full explicit `BDF1` assembly now returns remotely.
- Warm explicit full-physics assembly took `0.562 s`.
- Warm autodiff assembly exceeded `180 s`.
- Explicit `BDF1` speedup exceeds `320x` lower bound.
- Next focus: solver-level benchmarking, not Jacobian algebra.
