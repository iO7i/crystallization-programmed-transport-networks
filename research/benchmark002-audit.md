# Benchmark 002 evidence audit

Date: 2026-09-10

## Reproduction

- Canonical source revision: `32757a31d697025f7f84c06b7ce5a793810d96a1`.
- Committed result/evidence revision: `e10993ae28c4d85f968874baf867efea64829517`.
- Fresh clone: the public `benchmark-002` branch at the evidence revision, with a newly
  created environment.
- Canonical environment: Python 3.12.12, NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.1.
- Fresh-clone environment: Python 3.12.10, NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.1.
- Exact command: `python -m crystal_transport.cli benchmark run benchmark_002 --resolution 16 --output results/benchmark_002`.
- Numeric comparison: relative tolerance `1e-10`, absolute tolerance `1e-12`; all
  machine-readable metrics passed.
- Bitwise comparison: all CSV/JSON evidence files except the provenance manifest,
  all three figures, and the saved placement masks matched. The environment Python
  build string, manifest runtime, and public-clone source revision are expected
  provenance/environment differences; no bitwise claim is made for them.
- CI smoke command: `crystal-transport benchmark run benchmark_002 --resolution 8
  --validation-profile smoke --output ci-results/benchmark_002`. Its explicit smoke
  profile checks artifact generation, deterministic placement, no-op equivalence,
  denominator accounting, valid masks, and fail-closed inputs. It does not replace
  canonical numerical acceptance.

## Scientific checks

- Benchmark 001 zero-fraction baseline reproduced with maximum relative error `0.0`.
- `D_crystal = D_mobile` is an exact no-op control.
- Same seed reproduces the same placement digest.
- Maximum flux-balance error is below the declared `1e-5` tolerance.
- 002A uses three sweep seeds; 002B uses five placement seeds.
- 002B placement effects remain qualitatively visible across 16³, 24³, and 32³,
  while exact normalized values remain resolution-dependent.
- Matrix contrast `D_matrix = 1e-2, 1e-3, 1e-4` does not remove the placement effect.
- No connectivity transition is claimed: 002A and 002B preserve geometry by design.
- The earlier `rtol=1e-8` flux-balance failure remains recorded in
  `research/negative-results.md`; it is not presented as first-pass success.

## Interpretation boundary

Benchmark 002 establishes only the measured consequences of controlled, fixed-geometry
phenomenological crystalline-like mobility labels and placement rules in synthetic
voxel networks. The no-help result is a passive linear-diffusion control/property, not
a claim about all real crystallization. It does not establish real polymer crystallization, PES-b-PEO
performance, crystallization kinetics, crystal orientation physics, SCFT morphology,
inverse design, or experimental material performance.

The result did not demonstrate a beneficial mobility-only effect. Structural narrowing,
blockage, and any possible stabilization mechanism remain deferred to 002C, where the
help-versus-hurt question becomes open. The `baseline-flux-ranked` rule is retained only
as a literal operational hypothesis because it was not the most transport-degrading
placement in this run.
