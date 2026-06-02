# Run Metadata

Generated: 2026-06-02

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
