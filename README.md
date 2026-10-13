# LOGOS theory general · v0.5.0

Creator: **Vangelis D. Karalis**.

The associated manuscript has not yet been submitted. Cite the software using
[`CITATION.cff`](CITATION.cff); no journal publication or DOI is claimed.

Reproducible, dimensionless examples of whole-regimen design: dose amounts,
administration timing, loading–maintenance–coast input, periodic recovery and
repair, injury after cessation, recovery-gated constant-family comparisons and fixed-target minimum-dose frontiers. LOGOS means Longitudinal Optimization of
Global Outcomes and Strategies.

## Install and reproduce

The supported interpreter is **CPython 3.12 only**. Package metadata, local setup
and continuous integration use the same range. From an extracted checkout:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-lock.txt
python -m pip install --no-build-isolation --no-deps .
python -m unittest discover -s tests -v
logos-pkpd reproduce --output results
python scripts/check_reproducibility.py results
python scripts/verify_manifest.py
```

On Windows activate with `.venv\Scripts\activate`. Choose a new output directory
for every run. Results contain JSON verification, CSV tables/trajectories and
PNG/SVG figures. A failed numerical check terminates the command. The full
pipeline recomputes the four periodic cases, 27 local combinations and 18 stress
cases. It also checks selected periodic solutions by independent ODE integration,
finite benefit by integration across jumps, and 13 analytical benchmark checks.
Additional runs compute a full-constant-family recovery bound, a common-initial-state
comparison, deterministic parameter stress tests, imperfect-information values and
a four-target minimum-dose frontier. Companion analytical checks use distinct model assumptions and retain separate
names and results; they do not supply a clinical comparison.

The declared library versions are NumPy 2.3.5, SciPy 1.17.1 and Matplotlib 3.10.9.
A complete pinned installation is in `requirements-lock.txt`. Figures use the
bundled DejaVu Sans font for portable rendering. Numerical tolerances, rather
than byte identity of generated plots, determine reproducibility across systems.

## Synthetic data

`data/synthetic_pk.csv` contains **600 simulated measurements from 24 synthetic
subjects**, with dimensionless times, concentrations, elimination rates and
volumes. It contains no patient observations or fitted clinical profiles. The
seed is 20261003. Regenerate it with:

```bash
logos-pkpd generate-data --output data/synthetic_pk.csv
```

The illustrative cohort is separate from the normalized benchmark models;
it demonstrates a shareable input format and is not used to validate clinical
prediction, patient outcomes or a dosing policy. See `data/README.md`.

## Models and evidence

`docs/MODELS.md` states the equations, parameters, comparison contracts and
limitations. Finite values are cumulative benefit; periodic values are cycle
means. Every periodic candidate uses its own stationary initialization.
Periodic safety does not establish startup or cessation safety.

The implementation is an independent reconstruction of declared mathematical
models. Optimization provides numerical candidates within the stated families.
Only the ideal finite Emax loading–maintenance–coast solution has the stated
analytical global certificate. Sampled periodic searches are not global optimality
or interval-arithmetic certificates. The full mathematical injury tail is bounded
for the normalized finite example. These dimensionless examples are research
software and do not define drug-specific treatment recommendations.

## Repository layout

- `src/logos_pkpd/`: models, synthetic generator, CLI and independent audits.
- `tests/`: conservation, ODE comparison, negative-control and reproducibility tests.
- `scripts/`: release integrity and reference-result checks.
- `data/`: synthetic CSV and data dictionary.
- `reference/`: freshly generated release outputs and verification summary.
- `.github/workflows/`: clean CPython 3.12 installation and full reproduction.

## License and release

The MIT license applies to **all supplied code, documentation and synthetic
data**, permitting reproduction, modification and reuse subject to its notice.
Software version: **0.5.0**. See `CHANGELOG.md` and `docs/RELEASING.md`.

## Communication

This repository distributes research materials. Issues, discussions and pull
requests are disabled. Contributions and comments are not solicited.

GitHub topics:

`logos`, `pharmacokinetics`, `pharmacodynamics`, `pkpd`, `dose-optimization`,
`dosing-regimens`, `mathematical-modeling`, `computational-pharmacology`,
`model-informed-precision-dosing`, `treatment-optimization`, `drug-safety`.
