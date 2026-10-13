# Synthetic data dictionary

All quantities are dimensionless; identifiers denote simulated subjects.

| Column | Meaning |
| --- | --- |
| synthetic_id | SYN-001 through SYN-024 |
| time | 25 equally spaced observations on [0, 6] |
| dose | Single intravenous amount at time zero, 1 |
| k | Elimination rate, lognormal(log mean 0, log SD 0.15) |
| volume | Volume, lognormal(log mean 0, log SD 0.10) |
| true_concentration | dose / volume × exp(-k × time) |
| simulated_measurement | true concentration + N(0, 0.03²), clipped at zero |
| is_synthetic | Always true |

Generation: NumPy default_rng with seed 20261003, subject parameters then
measurement noise in ascending subject/time order. Code, CSV and this dictionary
are licensed under the repository MIT license. No real observations are included.
