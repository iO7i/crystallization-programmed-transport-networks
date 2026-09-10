# Benchmark 002 evidence audit

Date: 2026-09-10

## Reproduction

- Canonical source revision: `90ceafc640b5670c23eac846aad7f4786f6895f4`.
- Committed result/evidence revision: `8e456e6e7a09faaf25dffb14591f815f17257eeb`.
- Fresh clone: the public `benchmark-002` branch at the evidence revision.
- Environment: Python 3.12.12, NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.1.
- Exact command: `python -m crystal_transport.cli benchmark run benchmark_002 --resolution 16 --output results/benchmark_002`.
- Numeric comparison: relative tolerance `1e-10`, absolute tolerance `1e-12`; all
  machine-readable metrics passed.
- Bitwise comparison: all CSV/JSON evidence files except the provenance manifest,
  all three figures, and the saved placement masks matched. The manifest's runtime and
  public-clone source revision are expected provenance differences.

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

## Interpretation boundary

Benchmark 002 establishes only the measured consequences of controlled,
phenomenological crystalline-like mobility labels and placement rules in synthetic
voxel networks. It does not establish real polymer crystallization, PES-b-PEO
performance, crystallization kinetics, crystal orientation physics, SCFT morphology,
inverse design, or experimental material performance.

The result did not demonstrate a beneficial mobility-only effect. Structural narrowing,
blockage, and any possible stabilization mechanism remain deferred to a later bounded
experiment. The flux-ranked backbone rule is retained as a hypothesis because it was
not the most transport-degrading placement in this run.
