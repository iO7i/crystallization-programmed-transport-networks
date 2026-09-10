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
