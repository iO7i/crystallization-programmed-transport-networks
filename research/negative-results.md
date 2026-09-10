# Negative results and audit corrections

## Grid-sensitive exploratory fixture

The periodic Gaussian-random-field fixture does not preserve its connected-component
count across the tested resolutions: the phase has 5, 1, and 2 components at 16³,
24³, and 32³, respectively. It percolates in all three directions at all tested
resolutions, so the bounded Benchmark 001 conclusion survives without using this
fixture as a headline quantitative ranking. The fixture remains in the repository as
an explicitly exploratory, grid-sensitive result.

## Engineering corrections

An analytically calculated cylinder radius did not preserve the requested voxelized
phase fraction at coarse resolution. The generator now uses a deterministic rank
threshold and reports the measured fraction. This is a discretization correction, not
evidence about polymer physics.

An earlier transport run using solver `rtol=1e-8` failed the declared flux-balance
criterion in part of the benchmark. The benchmark solver was tightened to `rtol=1e-10`
in commit `d903f60`; the final manifest records both the solver tolerance and the
`1e-5` flux-balance acceptance criterion. The earlier failure is retained here rather
than presented as first-pass solver success.

## Benchmark 002 negative and unexpected results

In the mobility-only 002A model, every tested `D_crystal / D_mobile <= 1` condition
reduced or preserved transport; no beneficial transport result was found. This is the
expected passive linear-diffusion control/property of a fixed geometry with pointwise
positive mobility reduced or unchanged. It does not demonstrate that real crystallization
always hurts transport. Morphology/connectivity changes are deferred to 002C, where the
help-versus-hurt question becomes genuinely open.

The 002B placement result was not the initially intuitive ranking: on the gyroid,
seeded random placement produced the largest mean transport loss at the declared 20%
fraction and 0.1 mobility ratio, while the `baseline-flux-ranked` placement was less
damaging. The baseline-flux ranking therefore remains an explicit operational placement
hypothesis, not a validated graph or universal criticality measure.

No connectivity transition occurred because 002A and 002B intentionally preserve
geometry. Structural narrowing/blockage was deferred rather than inventing a positive
"stabilization" effect.
