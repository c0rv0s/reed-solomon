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

## Explicit Smooth-Domain Lower Bound

Let `H <= F_p^*` be a multiplicative subgroup of size `n`, let `M <= H` have size `ell`, and let
`S` be a union of `r` cosets of `M`, so `s = r * ell`. The vanishing polynomial

```text
P_S(x) = prod_{alpha in S} (x - alpha)
```

has only powers congruent to `0 mod ell`, because each coset contributes a factor of the form
`x^ell - y_i`. Therefore `x^s - P_S(x)` has degree at most `s - ell`.

For any received polynomial `R(x) = x^s + h(x)` with `deg h < k`, the polynomial

```text
f_S(x) = R(x) - P_S(x)
```

has degree `< k` whenever `s - ell < k`, and it agrees with `R` on every point of `S`.
Consequently,

```text
|Lambda(RS[H,k], (n-s)/n, R)| >= binom(n/ell, s/ell)
```

whenever `ell | n`, `ell | s`, `s > k`, and `s - ell < k`.

This exactly explains the validated triage row

```text
p = 193, n = 12, k = 3, s = 4, ell = 2, R(x) = x^4 + x^2
```

which gives `binom(6,2) = 15` nearby codewords. The exact subset count is also `15`, so this row is
a construction rather than a sampling artifact.

The Stage 6 grid in `reports/lower-bound-grid.csv` checks the requested power-of-two slice
`q_bits in {64,128,192,256}`, `eps_bits = 128`, `n in {2^10,2^12,2^16,2^20}`, and
`rho in {1/2,1/4,1/8,1/16}`. For `q_bits = 256`, no row beats the `eps*q` budget; the best margin
is `log2(list) - log2(eps*q) = 3 - 128 = -125`. In this parameter slice, the construction is a
real smooth-domain obstruction to account for, but not yet a prize-budget counterexample.

## Missing Prize-Level Lemma

Candidate upper bound:

```text
If H_{Q_m}(delta) <= 1 - R_m + log_{Q_m}(eps * q) / N_m - eta,
then every smooth-domain Reed-Solomon C_m satisfies
|Lambda(C_m, delta)| <= eps * q,
with explicit finite-size slack depending on eta, n, q, m, and rho.
```

The volume bound supplies the converse. This missing lemma is the proof target unless experiments
find smooth-domain counterexamples below the capacity line. The upper bound should be stated above
all explicit lower-bound families, including the coset-union construction above.
