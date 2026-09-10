# Scientific software audit (2026-09-10)

| Candidate | What it solves | Representation / applicability | License / maintenance signal | M1 decision |
|---|---|---|---|---|
| [PoreSpy](https://github.com/PMEAL/porespy) | 3-D image generation, characterization, pore-network extraction, and porous-media simulations | Strong voxel/image tooling; useful later for cross-checks and imported morphologies, but not a crystallization model | MIT; repository has tests, docs, and a current lockfile | Do not depend on it yet; keep direct NumPy/SciPy fixtures auditable |
| [OpenPNM](https://github.com/PMEAL/OpenPNM) | Pore-network modeling and transport on extracted networks | Network representation is valuable for later adapters, but extraction changes the object being solved | MIT; active public repository and documented install/citation path | Defer until a voxel-to-network comparison is a research question |
| [PSCF++](https://github.com/dmorse/pscfpp) | SCFT and partial-saddle-point field-theoretic simulation for periodic block-polymer systems | Highly relevant to later morphology generation; does not by itself provide crystallization kinetics | GPL; source build, GSL/FFTW, optional CUDA; maintained project history | Do not add as an M1 dependency; evaluate at M5 |
| [langevin-fts](https://github.com/yongdd/langevin-fts) | Langevin field-theoretic simulation | Potentially relevant to fluctuations and dynamics; outside the validated M1 scope | Audit deferred pending a reproducible current build and license check | Not used |
| [scft-id](https://github.com/olsenlabmit/scft-id) | SCFT and inverse design of block-polymer sequences/structures | Relevant only after forward morphology and transport are validated | MIT; research repository with environment and data assumptions | Inverse design is explicitly locked out of M1 |
| [FiPy](https://pages.nist.gov/fipy/en/latest/) | Finite-volume PDE solver with diffusion and phase-field capabilities | Could support a later phase-field model; its generality is unnecessary for the first transparent solver | NIST-maintained documentation and open source distribution | Not used in M1; direct sparse assembly makes the benchmark easier to audit |
| [freud](https://freud.readthedocs.io/en/latest/) | Particle/trajectory structure and local-environment analysis with periodic boxes | Useful for particle simulations, not for binary voxel transport in M1 | BSD-family distribution; current documentation exists | Defer |
| [MDAnalysis](https://docs.mdanalysis.org/) | Molecular-dynamics trajectory analysis | Relevant for future molecular transport validation, not for current synthetic fields | LGPLv3+; current docs and release information | Defer |

## Selection rationale

M1 uses NumPy, SciPy, and Matplotlib only. The morphology and transport representations
are small enough to inspect directly, and the solver is validated against homogeneous
and laminate limits before any external morphology or pore-network conversion is
introduced. Reuse remains possible, but package output will not be treated as validated
evidence without a repository-level comparison.

