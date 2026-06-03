# Run Metadata

Generated: 2026-06-03

## Environment

- Python: `3.13.3`
- SageMath: not installed (`sage` and `sageMath` were not found on `PATH`)
- pytest: not installed, so tests were run with stdlib `unittest`

## Commands

```text
env PYTHONPATH=src python3 -m unittest discover -s tests -v
env PYTHONPATH=src python3 -m rs_grand_list_decoding.rs_capacity_threshold
env PYTHONPATH=src python3 -m rs_grand_list_decoding.rs_capacity_threshold --csv reports/capacity-table.csv --summary reports/capacity-summary.md
env PYTHONPATH=src python3 -m rs_grand_list_decoding.exact_search --csv reports/exact-sweep.csv --summary reports/exact-sweep-summary.md
env PYTHONPATH=src python3 -m rs_grand_list_decoding.structured_bad_centers --csv reports/structured-centers.csv --summary reports/structured-centers-summary.md --comparison-csv reports/domain-comparison.csv --comparison-summary reports/domain-comparison-summary.md
env PYTHONPATH=src python3 -m rs_grand_list_decoding.triage_search --csv reports/triage-search.csv --summary reports/triage-search-summary.md --comparison-csv reports/triage-domain-comparison.csv --comparison-summary reports/triage-domain-comparison-summary.md --patterns-csv reports/triage-patterns.csv --patterns-summary reports/triage-patterns-summary.md --patterns-smooth-coset-csv reports/triage-patterns-smooth-coset.csv --patterns-smooth-coset-summary reports/triage-patterns-smooth-coset-summary.md --patterns-random-csv reports/triage-patterns-random.csv --patterns-random-summary reports/triage-patterns-random-summary.md
env PYTHONPATH=src python3 -m rs_grand_list_decoding.triage_validation --input reports/triage-search.csv --csv reports/triage-validation.csv --summary reports/triage-validation-summary.md
env PYTHONPATH=src python3 -m rs_grand_list_decoding.lower_bound_constructions --csv reports/lower-bound-grid.csv --summary reports/lower-bound-grid-summary.md
```

## Outputs

- `reports/test-results.txt`: saved unit-test output
- `reports/capacity-table.csv`: full capacity/volume-converse table with 512 data rows
- `reports/capacity-summary.md`: headline folded-mode table for `q = 2^256`, `n = 2^20`
- `reports/exact-sweep.csv`: first tiny-parameter exact max-list sweep
- `reports/exact-sweep-summary.md`: exact-sweep summary with successful/refused row counts
- `reports/structured-centers.csv`: first monomial/binomial structured-center sweep
- `reports/structured-centers-summary.md`: top structured-center counts
- `reports/domain-comparison.csv`: smooth-vs-random best-count comparison
- `reports/domain-comparison-summary.md`: top smooth-over-random ratios
- `reports/triage-search.csv`: sampled larger-field scale-and-triage search
- `reports/triage-search-summary.md`: top sampled triage rows
- `reports/triage-domain-comparison.csv`: triage smooth-vs-random comparison
- `reports/triage-domain-comparison-summary.md`: top triage smooth-over-random ratios
- `reports/triage-patterns.csv`: modular pattern groups from high triage rows
- `reports/triage-patterns-summary.md`: top recurring modular patterns
- `reports/triage-patterns-smooth-coset.csv`: subgroup/coset modular pattern groups
- `reports/triage-patterns-smooth-coset-summary.md`: top subgroup/coset modular patterns
- `reports/triage-patterns-random.csv`: random-domain raw and modulo-`p-1` pattern groups
- `reports/triage-patterns-random-summary.md`: top random-domain patterns
- `reports/triage-validation.csv`: validation reruns for top smooth triage rows
- `reports/triage-validation-summary.md`: empirical p-value validation summary
- `reports/lower-bound-grid.csv`: numerical scan of explicit coset-union lower bounds
- `reports/lower-bound-grid-summary.md`: budget/capacity summary for the lower-bound grid
