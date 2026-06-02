# Proof Ledger

## Convention Questions

- What exactly is `C^{equiv m}` in the prize statement: folded RS, direct interleaving, or another
  interleaved convention?
- Is the intended regime fixed `n` with `q` sufficiently large, or simultaneous growth of `n` and
  `q` under SNARK-relevant constraints?
- Are balls closed (`<= delta n`) or strict (`< delta n`) at boundary radii?

## Known Converse

For any code `C_m` with block length `N_m`, alphabet size `Q_m`, and rate `R_m`, averaging over
received words gives:

```text
|C_m| * V_{Q_m}(N_m, t) > eps * q * Q_m^{N_m}
    => max_r |Lambda(C_m, t / N_m, r)| > eps * q.
```

This volume bound is code-independent and is therefore a rigorous upper limit on any proposed
`delta_C*`.

## Missing Prize-Level Lemma

Candidate upper bound:

```text
If H_{Q_m}(delta) <= 1 - R_m + log_{Q_m}(eps * q) / N_m - eta,
then every smooth-domain Reed-Solomon C_m satisfies
|Lambda(C_m, delta)| <= eps * q,
with explicit finite-size slack depending on eta, n, q, m, and rho.
```

The volume bound supplies the converse. This missing lemma is the proof target unless experiments
find smooth-domain counterexamples below the capacity line.

