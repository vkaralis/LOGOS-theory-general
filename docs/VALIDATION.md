# Release validation

Locally executed in a newly created virtual environment, without system site
packages, using CPython 3.12.14 on macOS arm64. Dependencies were downloaded and
installed from requirements-lock.txt. The package was built and installed with
`pip install --no-build-isolation --no-deps .`; `pip check` found no broken requirements.
No private source archive, workbook or document is needed at runtime.

- Six mathematical/data tests passed.
- All 24 seeded timing optimization starts succeeded in the final source.
- 45 periodic reference checks passed (four cases, 27 local and 18 stress cells).
- 13 independent analytical benchmark checks passed.
- 73 recovery comparison, parameter stress, information and finite-bound checks passed.
- 20 minimum-dose frontier checks passed.
- Independent finite ODE maximum objective error: 6.61e-13 or less.
- Selected periodic ODE mean-metric maximum error: 4.36e-13 or less.
- A separate complete run passed 3052 reference comparisons and exact synthetic
  CSV identity. PNG/SVG figures were opened for visual inspection.

The release reference is a fresh computation from the packaged source.
Reproduction uses numerical tolerances; numerical candidates and floating-point
bounds do not establish additional global or clinical guarantees.

Local execution used macOS arm64. The complete pipeline also passed on GitHub
Actions in a fresh Linux CPython 3.12 environment on 2026-10-13:
[verified run](https://github.com/vkaralis/LOGOS-theory-general/actions/runs/38446883227).
The run passed file-integrity verification, all six tests, all scientific audits
and 3,052 reference comparisons. Other Python minor versions are unsupported.

The release reproduces the dimensionless models and synthetic cohort described
in the public documentation. Clinical ANC prediction, fitted patient-specific
warfarin/remifentanil comparisons and historical implementation statistics are
outside this release's validated scope.
