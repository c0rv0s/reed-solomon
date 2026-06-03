# Triage Validation Summary

Reruns top smooth-domain triage candidates at larger subset budgets, then compares the
same center family against deterministic random domains with matching parameters.

- Input: `reports/triage-search.csv`
- Rows: `10`
- Smooth budgets: `100, 1000, 10000`

| p-value | smooth best | random max | random mean | p | n | k | s | type | parameters | low-k explained |
|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| 0.000 | 15 | 6 | 2.490 | 193 | 12 | 3 | 4 | binomial | `a=2;b=4;lambda=1` | False |
| 0.000 | 6 | 2 | 0.420 | 193 | 12 | 1 | 2 | monomial | `exponent=2` | True |
| 0.000 | 6 | 2 | 0.420 | 193 | 12 | 1 | 2 | binomial | `a=0;b=2;lambda=1` | False |
| 0.000 | 6 | 2 | 0.420 | 193 | 12 | 1 | 2 | shifted_binomial | `c=1;a=0;b=2;lambda=1` | False |
| 0.000 | 6 | 3 | 0.800 | 97 | 12 | 1 | 2 | monomial | `exponent=2` | True |
| 0.000 | 6 | 3 | 0.800 | 97 | 12 | 1 | 2 | binomial | `a=0;b=2;lambda=1` | False |
| 0.000 | 6 | 3 | 0.800 | 97 | 12 | 1 | 2 | shifted_binomial | `c=1;a=0;b=2;lambda=1` | False |
| 0.000 | 5 | 2 | 0.290 | 193 | 8 | 2 | 3 | trinomial | `a=2;b=3;c=4;lambda1=1;lambda2=1` | False |
| 0.000 | 5 | 2 | 0.710 | 193 | 12 | 1 | 2 | binomial | `a=2;b=4;lambda=1` | False |
| 0.000 | 5 | 2 | 0.710 | 193 | 12 | 1 | 2 | shifted_binomial | `c=1;a=2;b=4;lambda=1` | False |
