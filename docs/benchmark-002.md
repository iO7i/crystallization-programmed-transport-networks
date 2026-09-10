# Benchmark 002 — When Crystallization Helps or Hurts Transport

Benchmark 002 is a phenomenological crystalline-like mobility-perturbation benchmark. It starts
from the stable Benchmark 001 fixtures and labels selected voxels of the existing
transport phase as crystalline-like low-mobility subdomains. It does not simulate
polymer crystallization, nucleation, growth, chain folding, SCFT, molecular dynamics,
or a specific polymer system.

002A and 002B are fixed-geometry passive linear-diffusion experiments. They change only
the positive pointwise diffusivity field `D(x)` while leaving morphology, boundary
conditions, and transport equations fixed. Therefore a mobility-only increase above the
unperturbed effective transport is not an expected mechanism in this model. The observed
no-help result is a validated control/property of this model, not evidence that real
crystallization can never improve transport. A genuine help-versus-hurt question begins
in 002C, where morphology/connectivity or another physical term may be changed.

## States and denominator

Each voxel has one of three transport states:

```text
matrix                 D_matrix = 1e-3
mobile transport       D_mobile = 1
crystalline-like       D_crystal = ratio * D_mobile
```

The crystalline-like fraction is always defined as:

```text
selected crystalline-like voxels / existing transport-phase voxels
```

The output also records the selected count, actual transport-phase fraction, and total
material fraction. Geometry is unchanged in 002A and 002B.

## Experiments

002A sweeps crystalline-like fractions `0, 0.1, 0.2, 0.4, 0.6` and mobility ratios
`1, 0.1, 0.01, 0.001` using seeded random placement on lamellae, cylinders, and the
analytical gyroid. Three seeds are used for the sweep.

002B fixes the fraction at 20% of the transport phase and the mobility ratio at 0.1
on the gyroid, comparing five seeds across four algorithmic placement modes:

- `random`: seeded uniform sampling without replacement;
- `interface`: smallest periodic Euclidean distance to the matrix interface;
- `baseline-flux-ranked`: largest baseline local flux-magnitude score; this is a literal
  ranking rule, not a validated graph or topological criticality measure;
- `low_criticality`: smallest baseline local flux-magnitude score.

The baseline flux score combines the three unperturbed directional solutions. It is
computed before the perturbation and is never defined from the perturbed solution.

Structural narrowing/blockage experiments are intentionally deferred. Mobility-only
perturbations cannot create a connectivity transition, but connectivity is recorded
for every condition to verify that this remains true.

## Controls and validation

- The zero-fraction control is compared directly with the committed Benchmark 001
  directional result at 32³.
- `D_crystal = D_mobile` is a no-op transport control.
- Repeated placement with the same seed must have the same selection digest.
- A crystalline-like label with equal mobility must leave the voxelwise field unchanged.
- Independent low-face/high-face fluxes must satisfy the existing `1e-5` tolerance.
- Representative baseline, mobility, and placement cases are checked at 16³, 24³,
  and 32³.

## Outputs

`results/benchmark_002/` contains `runs.csv`, aggregate summaries, grid and contrast
tables, controls, manifest/environment records, and figures. Figures are regenerated
from the written result tables and the saved placement masks.

The CLI defaults to the `canonical` validation profile. CI uses the explicit
`--validation-profile smoke` profile at `8³`; that profile checks artifact generation,
deterministic placement, no-op equivalence, phase-denominator accounting, valid masks,
and fail-closed inputs. It does not replace canonical numerical acceptance, which is
run on the declared Benchmark 002 grids.

The Gaussian random field is excluded from headline Benchmark 002 conclusions because
its connected-component count was not resolution-stable in Benchmark 001. The gyroid
remains a first-harmonic analytical level-set phase, not a double gyroid or evidence of
physical block-copolymer self-assembly.
