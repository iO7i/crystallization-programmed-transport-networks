# Benchmark 001

## Question

If the measured conducting-phase fraction and intrinsic phase coefficients are held
fixed, how much can morphology, connectivity, and orientation change directional
effective diffusivity?

## Model

The benchmark solves the steady scalar diffusion equation

```text
∇ · (D(x) ∇c) = 0
```

on a cell-centred Cartesian voxel grid. The transport phase has `D_transport = 1.0`;
the default matrix has `D_matrix = 1e-3`. This is a near-insulating model assumption,
not an assertion that a polymer matrix has zero mobility.

For each measured axis, `c = 1` and `c = 0` are applied at the two external faces;
the other two axes are periodic. Harmonic face conductances are used at phase
interfaces. The reported scalar is the area-averaged boundary flux multiplied by
the sample length and divided by the imposed concentration difference.

The implementation is intentionally independent of SCFT and crystallization models.
The four M1 fixtures are:

- `lamellae`: planar binary layers, with the layer normal along x;
- `cylinders`: a square lattice of four parallel conducting cylinders, with axes z;
- `gyroid_level_set`: a thresholded first-harmonic mathematical gyroid reference;
- `gaussian_random_field`: a periodic low-frequency Gaussian random-field reference.

The last two names describe the mathematical construction only. Neither is presented
as an SCFT-predicted or crystallized polymer morphology.

All fixtures use deterministic rank thresholding to select the requested number of
voxels from their generator field. The result records the rank field, threshold, and
selected voxel count. This controls the measured fraction without claiming that the
thresholded geometry is equivalent to a continuum morphology. At coarse resolutions,
a fractional voxel plane can change a nominal laminate geometry, so the benchmark
reports geometry and topology at every resolution rather than hiding that effect.

The gyroid is exactly one thresholded first-harmonic level-set phase in a periodic
mathematical field. It is not labelled as single gyroid, double gyroid, or an SCFT
solution. It does not establish physically self-assembled block-copolymer gyroid
behavior.

## Validation

Before the complex fixtures, the same solver is checked against:

1. a uniform medium in x, y, and z;
2. the series effective coefficient of a 50/50 laminate normal to the layers;
3. the parallel effective coefficient of the same laminate.

The validation results are written to `results/benchmark_001/validation.json`.

Independent low-face and high-face fluxes are compared. The acceptance tolerance is
`1e-5` relative flux difference and is recorded in `manifest.json` and each transport
row. The solver residual is also retained.

## Executed result

The definitive clean 32³ result is in `results/benchmark_001/transport_metrics.csv`.
The file reports `D_eff,x`, `D_eff,y`, and `D_eff,z` as separate directional rows,
along with fluxes and flux-balance errors. `morphology_metrics.csv` records the target
and actual fractions, generator thresholds, connected components, percolation, area
density, and mean phase chord lengths.

The benchmark’s central result is not a morphology ranking. It is that matched phase
fraction does not determine directional transport under the declared scalar-diffusion
model. The generated figure is a compact view; CSV and JSON files are authoritative.

## Resolution sensitivity

The committed study compares `16³`, `24³`, and `32³`. The qualitative directional
states of lamellae and cylinders remain stable; the gyroid remains one component and
percolates in all three axes. The Gaussian random-field fixture percolates in all
tested axes, but its component count is not stable across resolution. Interfacial
area density, voxel chord lengths, and some absolute transport values also move with
resolution. This is a benchmark-design result, not a reason to tune the generator for
a cleaner ranking narrative.

The detailed per-resolution geometry/topology record is
`morphology_convergence.csv`; the interpreted summary is `stability_summary.json`.

## Transport-contrast sensitivity

At `16³`, the benchmark sweeps `D_matrix = 10^-2, 10^-3, 10^-4` with
`D_transport = 1`. Results are in `contrast_sensitivity.csv`. This keeps the matrix
coefficient explicit and tests whether conclusions depend on the near-insulating
contrast assumption.

## Reproduction

From a clean environment:

```bash
python -m pip install -e ".[dev]"
crystal-transport validate
crystal-transport benchmark run benchmark_001 --resolution 32 --output results/benchmark_001
```

The command regenerates CSV metrics, JSON provenance, and figures directly from the
source code. The result manifest records the configuration digest, source revision,
environment, solver and flux tolerances, input coefficients, and whether tracked
source files were clean before execution.

## Limitations

- synthetic voxel geometries are not experimental polymer morphologies;
- no crystallization, chain architecture, segmental mobility, chemistry, or temperature
  dependence is modeled;
- scalar diffusion is not ionic conductivity, permeability, solvent transport, or gas
  transport;
- finite resolution changes interfaces and can change metrics;
- the matrix coefficient is a declared model parameter, not calibrated material data;
- the periodic percolation metric is a face-to-face component test;
- no experimental calibration or validation is claimed.

