# Crystallization-Programmed Transport Networks

Computational research into how polymer morphology, topology, confinement, and eventually crystallization history may determine transport networks.

```text
polymer architecture / processing
             ↓
       3-D morphology
             ↓
      network topology
             ↓
          transport

crystallization is a future coupling variable, not an M1 model
```

## Research question

Can molecular topology, spatial confinement, and crystallization history become design variables for predicting the connectivity, dimensions, and transport properties of self-assembled polymer networks?

The project does not assume that a universal design rule exists. It tests which relationships are reproducible, which are material-specific, and which are too uncertain for prediction.

## What has actually been done

Benchmark 001 is an executed morphology-to-transport experiment using four deterministic 3-D voxel fixtures at matched measured phase fraction. A validated finite-volume solver measures directional steady scalar diffusion in x, y, and z.

Reproduce it with:

```bash
python -m pip install -e ".[dev]"
crystal-transport benchmark run benchmark_001 --resolution 32 --output results/benchmark_001
```

The result figure is [Benchmark 001 overview](results/benchmark_001/figures/benchmark_001_overview.png), with source data in [transport_metrics.csv](results/benchmark_001/transport_metrics.csv) and provenance in [manifest.json](results/benchmark_001/manifest.json).

## Evidence status

- **REPRODUCED:** morphology and orientation change effective transport under Benchmark 001 controls.
- **REPRODUCED:** a clean-clone 32³ rerun regenerated all outputs; metric files and figures match the committed reference, with numeric comparison tolerance recorded in `results/benchmark_001/reproduction_audit.json`.
- **SUPPORTED:** confinement and crystallization history can alter morphology in published crystalline-polymer systems.
- **OPEN:** whether crystallization history can become a predictive design variable for transport topology.
- **NOT CLAIMED:** polymer crystallization, PES-b-PEO performance, real membrane performance, SCFT-generated morphology, inverse design, experimental validation, or universal transport rules.

See the [claim ledger](research/claim-ledger.csv), [verified lineage](research/lineage/ikehara.md), [frontier note](research/frontier.md), and [software audit](research/software-audit.md).

## Benchmark 001 result

At `32³`, every fixture has the same measured transport-phase fraction to voxel-count precision. The directional result table, resolution study, contrast study, component/percolation metrics, chord lengths, flux-balance errors, and generator thresholds are generated under `results/benchmark_001/`.

The benchmark’s central result is not a morphology ranking. It is that matched phase fraction does not determine directional transport under the declared scalar-diffusion model. The resolution study also shows that absolute values and some geometric descriptors remain grid-sensitive. The Gaussian-random-field component count is not stable across the tested resolutions and is therefore exploratory, not a headline quantitative result.

## Benchmark 002 — phenomenological crystalline-like perturbations

Benchmark 002 reuses the stable lamella, cylinder, and analytical gyroid fixtures. It
labels a controlled fraction of the existing transport phase as a crystalline-like
low-mobility subdomain and measures how mobility penalty and spatial placement affect
transport. It does not model real crystallization, PES-b-PEO behavior, SCFT, or
experimental material performance.

Run it with:

```bash
crystal-transport benchmark run benchmark_002 --resolution 16 --output results/benchmark_002
```

The denominator is always the existing transport-phase voxel count. The first release
scope contains mobility-only 002A and spatial-placement 002B; structural narrowing and
blockage are deferred until those experiments are stable. See
[docs/benchmark-002.md](docs/benchmark-002.md) and
[benchmarks/benchmark_002.yaml](benchmarks/benchmark_002.yaml).

## Validation and limitations

The solver passes uniform-medium and laminate series/parallel checks before running Benchmark 001. Independent low-face/high-face fluxes are compared with an acceptance tolerance recorded in the manifest. The matrix coefficient is explicit and is swept across `D_matrix = 10⁻², 10⁻³, 10⁻⁴` at `16³`; it is a near-insulating model assumption, not hidden physics.

The benchmark models scalar diffusion only. It has no polymer chemistry, chain mobility, crystallization kinetics, temperature dependence, ionic conductivity, permeability, or experimental data. Detailed assumptions and limitations are in [docs/benchmark-001.md](docs/benchmark-001.md).

The gyroid fixture is exactly a thresholded first-harmonic analytical gyroid level-set geometry. It represents one binary level-set phase in a periodic mathematical field; it is not a double-gyroid SCFT morphology and does not establish physically self-assembled block-copolymer behavior.

## Research lineage and independence

This project was initiated by Hosam Al-Khairat as a return to questions originating in earlier chemical-engineering research on polymer crystallization, now approached through scientific computing. The earlier unpublished measurements are not used as validated reference data.

This is an independent research-software project and does not represent Kanagawa University or the Ikehara Laboratory. No university or laboratory endorsement, collaboration, or affiliation is implied.

## Repository structure

```text
src/crystal_transport/       minimal morphology, connectivity, transport, and CLI code
tests/                       scientific-invariant and reproduction tests
benchmarks/                  machine-readable experiment definition
research/                    claim ledger, literature, lineage, and software audit
docs/                        benchmark assumptions and interpretation
results/benchmark_001/       compact evidence and generated figures
```

Future crystallization, SCFT, inverse-design, and experimental directories are intentionally absent until those models exist and are validated.

## Development checks

```bash
ruff check src tests
mypy
pytest
python -m build
```

The project uses Python 3.12 with NumPy, SciPy, and Matplotlib. A small CI smoke benchmark is defined in [.github/workflows/ci.yml](.github/workflows/ci.yml). The v0.1.0-alpha release candidate has passed the fresh-clone reproduction and final citation, privacy, and affiliation audit; the reproduced result records the source and evidence revisions explicitly.

## Research principle

This project does not assume that a universal crystallization-to-transport design rule exists. A negative result that establishes a boundary of predictability is a useful result.
