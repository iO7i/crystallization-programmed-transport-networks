# v0.1.0-alpha release audit

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

The Gaussian-random-field fixture remains in the result package as an explicitly
exploratory, grid-sensitive case. It is excluded from headline quantitative conclusions
until its geometry is stabilized or its role is otherwise resolved.
