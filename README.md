# RS Grand List-Decoding

Research harness for attacking the Reed-Solomon grand list-decoding challenge from the Proximity
Prize.

## Current Artifacts

- `SPEC.md`: project spec, mathematical target, and implementation plan.
- `docs/proof-ledger.md`: proof target and convention questions.
- `src/rs_grand_list_decoding/`: pure-Python threshold and prime-field utilities.
- `tests/`: stdlib `unittest` test suite.
- `reports/capacity-table.csv`: saved capacity/volume-converse experiment results.
- `reports/capacity-summary.md`: headline summary for `q = 2^256`, `n = 2^20`.
- `reports/exact-sweep.csv`: first tiny-parameter exact max-list sweep.
- `reports/exact-sweep-summary.md`: exact-sweep summary and refusal counts.
- `reports/structured-centers.csv`: first structured monomial/binomial bad-center sweep.
- `reports/structured-centers-summary.md`: top structured-center counts.
- `reports/domain-comparison.csv`: smooth-vs-random best-count comparison.
- `reports/domain-comparison-summary.md`: top smooth-over-random ratios.
- `reports/triage-search.csv`: sampled scale-and-triage search for larger prime fields.
- `reports/triage-domain-comparison.csv`: triage smooth-vs-random comparison.
- `reports/triage-patterns.csv`: modular exponent pattern groups from high triage rows.

## Run

```text
env PYTHONPATH=src python3 -m unittest discover -s tests -v
env PYTHONPATH=src python3 -m rs_grand_list_decoding.rs_capacity_threshold
env PYTHONPATH=src python3 -m rs_grand_list_decoding.rs_capacity_threshold --csv reports/capacity-table.csv --summary reports/capacity-summary.md
env PYTHONPATH=src python3 -m rs_grand_list_decoding.exact_search --csv reports/exact-sweep.csv --summary reports/exact-sweep-summary.md
env PYTHONPATH=src python3 -m rs_grand_list_decoding.structured_bad_centers --csv reports/structured-centers.csv --summary reports/structured-centers-summary.md --comparison-csv reports/domain-comparison.csv --comparison-summary reports/domain-comparison-summary.md
env PYTHONPATH=src python3 -m rs_grand_list_decoding.triage_search --csv reports/triage-search.csv --summary reports/triage-search-summary.md --comparison-csv reports/triage-domain-comparison.csv --comparison-summary reports/triage-domain-comparison-summary.md --patterns-csv reports/triage-patterns.csv --patterns-summary reports/triage-patterns-summary.md
```

## Notes

The baseline package is pure Python. SageMath can be added later as an optional backend for
extension fields, interpolation, and larger exact finite-field experiments.
