# v0.2.0-alpha release audit

Date: 2026-09-10

## Citation and provenance

- The release metadata is in `CITATION.cff` and names only the project author.
- The canonical Benchmark 001 result was generated from source revision
  `d903f60c1baa108f35dfcb51e1f66085fe364adc`.
- The committed reference result was added in evidence revision
  `00e5e4adaac97e03b9fed154ad2399d78ac1a5c8`; the reproduction record is
  `results/benchmark_001/reproduction_audit.json`.
- The reproduction used a fresh clone, a new Python 3.12 environment, the documented
  command, an empty result directory, and regenerated machine-readable outputs and
  figures. Metric files were numerically equivalent at relative tolerance `1e-10` and
  absolute tolerance `1e-12`; the metric CSVs, validation, environment file, and all
  figures were bitwise identical. `stability_summary.json` is semantically identical
  but key ordering is not bitwise identical; `manifest.json` necessarily changes its
  runtime field between runs.
- The canonical Benchmark 002 source revision is
  `32757a31d697025f7f84c06b7ce5a793810d96a1`; its committed result/evidence revision is
  `e10993ae28c4d85f968874baf867efea64829517`, with the final audit and claim-ledger
  update in the subsequent release-gate revision. The exact command, configuration
  digest, environment, solver tolerance, flux criterion, actual phase fractions,
  directional transport, grid study, contrast study, and limitations are recorded in
  `results/benchmark_002/manifest.json`, `validation.json`, and the accompanying CSVs.
- A fresh clone with a newly created environment ran the actual B002 command and the
  actual B001 32³ regression. B002 machine-readable metrics were numerically equivalent
  at relative tolerance `1e-10` and absolute tolerance `1e-12`; CSV/JSON evidence,
  figures, and placement masks were bitwise identical except for expected provenance
  and Python-build-string differences. The complete comparison is
  `results/benchmark_002/reproduction_audit.json`.

## Independence and affiliation

- This is an independent research-software project.
- The repository does not represent Kanagawa University or the Ikehara Laboratory and
  makes no endorsement, collaboration, or affiliation claim.
- Earlier unpublished measurements are not used as validated input data.
- Published lineage and citations are kept in `research/lineage/ikehara.md` and
  `research/bibliography/`; they are not presented as repository-generated results.

## Claim boundary

Benchmark 001 supports only the bounded claim that morphology and orientation can
change directional effective scalar diffusivity under the declared binary voxel model,
with the reported controls and numerical checks. It does not support claims about
polymer crystallization, PES-b-PEO material performance, real membrane performance,
SCFT-generated morphology, inverse design, or experimental validation.

Benchmark 002A/002B support only the tested fixed-geometry, passive linear-diffusion
property/control: pointwise mobility reductions or unchanged mobility did not improve
effective transport, and placement of the reduced-mobility labels changed transport.
They are phenomenological perturbation experiments, not a crystallization simulation;
they do not reconstruct morphology or establish whether real crystallization helps or
hurts transport. Structural help/hurt remains deferred to 002C. The
`baseline-flux-ranked` label is an operational ranking rule, not validated network
criticality.

The earlier `rtol=1e-8` flux-balance failure is retained in
`research/negative-results.md`; the release does not imply that every configuration
passed on the first attempt.

The Gaussian-random-field fixture remains in the result package as an explicitly
exploratory, grid-sensitive case. It is excluded from headline quantitative conclusions
until its geometry is stabilized or its role is otherwise resolved.
