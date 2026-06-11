# RS Grand List-Decoding

Research harness for attacking the Reed-Solomon grand list-decoding challenge from the Proximity
Prize.

## Current Artifacts

- `SPEC.md`: project spec, mathematical target, and implementation plan.
- `docs/proof-ledger.md`: proof target and convention questions.
- `docs/organizer-question.md`: short smooth-domain convention question to send to organizers.
- `docs/mixed-smooth-counterexample.md`: self-contained proof note for the mixed-smooth witness.
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
- `reports/triage-patterns-smooth-coset.csv`: subgroup/coset exponent patterns modulo `n`.
- `reports/triage-patterns-random.csv`: random-domain exponent patterns using raw exponents and modulo `p-1`.
- `reports/triage-validation.csv`: reruns top smooth triage candidates against 100 random domains.
- `reports/triage-validation-summary.md`: empirical p-value summary for validated candidates.
- `reports/lower-bound-grid.csv`: numerical scan of explicit smooth-domain coset-union lower bounds.
- `reports/lower-bound-grid-summary.md`: budget/capacity summary for the lower-bound grid.
- `reports/quotient-lower-bound.csv`: quotient-optimized mixed-smooth coset-union search.
- `reports/quotient-lower-bound-summary.md`: headline quotient search summary.
- `reports/quotient-materialized-primes.csv`: concrete prime fields for selected quotient candidates.
- `reports/mixed-smooth-counterexample.json`: machine-readable flagship counterexample witness.

## Run

```text
env PYTHONPATH=src python3 -m unittest discover -s tests -v
env PYTHONPATH=src python3 -m rs_grand_list_decoding.rs_capacity_threshold
env PYTHONPATH=src python3 -m rs_grand_list_decoding.rs_capacity_threshold --csv reports/capacity-table.csv --summary reports/capacity-summary.md
env PYTHONPATH=src python3 -m rs_grand_list_decoding.exact_search --csv reports/exact-sweep.csv --summary reports/exact-sweep-summary.md
env PYTHONPATH=src python3 -m rs_grand_list_decoding.structured_bad_centers --csv reports/structured-centers.csv --summary reports/structured-centers-summary.md --comparison-csv reports/domain-comparison.csv --comparison-summary reports/domain-comparison-summary.md
env PYTHONPATH=src python3 -m rs_grand_list_decoding.triage_search --csv reports/triage-search.csv --summary reports/triage-search-summary.md --comparison-csv reports/triage-domain-comparison.csv --comparison-summary reports/triage-domain-comparison-summary.md --patterns-csv reports/triage-patterns.csv --patterns-summary reports/triage-patterns-summary.md --patterns-smooth-coset-csv reports/triage-patterns-smooth-coset.csv --patterns-smooth-coset-summary reports/triage-patterns-smooth-coset-summary.md --patterns-random-csv reports/triage-patterns-random.csv --patterns-random-summary reports/triage-patterns-random-summary.md
env PYTHONPATH=src python3 -m rs_grand_list_decoding.triage_validation --input reports/triage-search.csv --csv reports/triage-validation.csv --summary reports/triage-validation-summary.md
env PYTHONPATH=src python3 -m rs_grand_list_decoding.lower_bound_constructions --csv reports/lower-bound-grid.csv --summary reports/lower-bound-grid-summary.md
env PYTHONPATH=src python3 -m rs_grand_list_decoding.quotient_lower_bound_search --csv reports/quotient-lower-bound.csv --summary reports/quotient-lower-bound-summary.md --materialized-csv reports/quotient-materialized-primes.csv --counterexample-json reports/mixed-smooth-counterexample.json
```

## Notes

The baseline package is pure Python. SageMath can be added later as an optional backend for
extension fields, interpolation, and larger exact finite-field experiments.
