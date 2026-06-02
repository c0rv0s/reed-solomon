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

## Run

```text
env PYTHONPATH=src python3 -m unittest discover -s tests -v
env PYTHONPATH=src python3 -m rs_grand_list_decoding.rs_capacity_threshold
env PYTHONPATH=src python3 -m rs_grand_list_decoding.rs_capacity_threshold --csv reports/capacity-table.csv --summary reports/capacity-summary.md
```

## Notes

The baseline package is pure Python. SageMath can be added later as an optional backend for
extension fields, interpolation, and larger exact finite-field experiments.

