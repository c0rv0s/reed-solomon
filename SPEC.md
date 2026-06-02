# Reed-Solomon Grand List-Decoding Spec

## Objective

Build a research/code workspace for attacking the Proximity Prize list-decoding challenge for
Reed-Solomon codes over smooth evaluation domains.

The primary goal is to test, refine, and either prove or falsify the capacity-threshold formula for
the grand list-decoding challenge:

For `C = RS[F, L, k]`, with `L` a smooth multiplicative subgroup or coset, rate
`rho = k / |L|` in `{1/2, 1/4, 1/8, 1/16}`, target `eps = 2^-128`, and constant `m`, determine the
largest radius `delta_C*` such that

```text
|Lambda(C^{equiv m}, delta_C*)| <= eps * |F|.
```

The codebase should support both computational experimentation and proof-focused bookkeeping.

## Research Thesis

The proposed target threshold is the `Q`-ary list-decoding capacity line, adjusted by the finite
`eps * |F|` list-size budget.

Let:

```text
q = |F|
n = |L|
k = rho * n
C_m = C^{equiv m}
N_m = block length of C_m
Q_m = alphabet size of C_m
R_m = log_{Q_m}(|C_m|) / N_m
```

Define the `Q`-ary Hamming ball volume:

```text
V_Q(N, t) = sum_{i=0}^t binom(N, i) * (Q - 1)^i
```

and the `Q`-ary entropy:

```text
H_Q(delta) =
    delta * log_Q(Q - 1)
  - delta * log_Q(delta)
  - (1 - delta) * log_Q(1 - delta).
```

The exact universal volume converse is:

```text
|C_m| * V_{Q_m}(N_m, t) > eps * q * Q_m^{N_m}
    => exists received word r with list size > eps * q at radius t / N_m.
```

Therefore no solution can exceed:

```text
delta_vol* =
    (1 / N_m) * max {
        t : V_{Q_m}(N_m, t) <= eps * q * Q_m^{N_m * (1 - R_m)}
    }.
```

The asymptotic capacity target is:

```text
H_{Q_m}(delta_cap*) =
    1 - R_m + log_{Q_m}(eps * q) / N_m
```

or:

```text
delta_cap* =
    H_{Q_m}^{-1}(1 - R_m + log_{Q_m}(eps * q) / N_m).
```

For folded mode, the expected convention is:

```text
Q_m = q^m
N_m = n / m
R_m = rho
```

so:

```text
delta_cap* =
    H_{q^m}^{-1}(1 - rho + (1 + log_q(eps)) / n).
```

With `eps = 2^-128`:

```text
delta_cap* =
    H_{q^m}^{-1}(1 - rho + (1 - 128 / log_2(q)) / n).
```

## Critical Ambiguity

The public prize statement does not fully define `C^{equiv m}` or the intended asymptotic regime.
This workspace must keep folded and direct interleaving separate.

Folded mode:

```text
N_m = n / m
Q_m = q^m
R_m = rho
```

Direct interleaving mode:

```text
N_m = n
Q_m = q^m
R_m = rho
```

This distinction changes the `m`-dependence of the entropy correction and may change which
formula is prize-relevant.

## Scope

This project should build tools and notes for the grand list-decoding challenge first. Mutual
correlated agreement is related, but should remain secondary unless a list-decoding route stalls.

In scope:

- Finite-field arithmetic for prime fields and extension fields when needed.
- Smooth multiplicative subgroup or coset construction.
- Scalar RS encoding.
- Folded RS encoding.
- Direct interleaved RS encoding.
- Exact brute-force list-size computation for tiny fields.
- Exact log-space volume thresholds for large symbolic parameters.
- Structured searches for bad received words over smooth domains.
- Experiment ledgers that connect observed failures/successes to candidate lemmas.

Out of scope for the first implementation pass:

- Efficient Guruswami-Sudan decoding as a core dependency.
- Large-scale exhaustive search beyond tiny toy fields.
- Proving MCA bounds.
- Building SNARK protocols.

## Modules

### `rs_capacity_threshold.py`

Pure Python numerical module.

Responsibilities:

- Compute `H_Q(delta)`.
- Invert `H_Q` by bisection.
- Compute log-space Hamming ball volumes.
- Compute exact finite volume-converse threshold.
- Compute entropy approximation threshold.
- Print comparison tables for `rho in {1/2, 1/4, 1/8, 1/16}`.

Required API:

```python
def h_q(delta: float, ln_q: float) -> float: ...
def inverse_h_q(target: float, ln_q: float, iters: int = 100) -> float: ...
def log_q_ball(N: int, t: int, ln_q: float) -> float: ...
def threshold_params(
    q_bits: float,
    n: int,
    rho: float,
    m: int,
    eps_bits: float = 128.0,
    mode: str = "folded",
) -> dict: ...
```

Acceptance:

- For `q = 2^256`, `n = 2^20`, `m = 1`, output capacity candidates approximately:
  - `rho = 1/2`: `0.49609`
  - `rho = 1/4`: `0.74683`
  - `rho = 1/8`: `0.87288`
  - `rho = 1/16`: `0.93618`
- For `m = 4`, outputs move closer to `1 - rho`.
- Exact grid threshold must never exceed the entropy threshold by more than one grid step plus
  numerical tolerance.

### `finite_field.py`

Small finite-field module for experiments.

Responsibilities:

- Implement arithmetic for `GF(p)` first.
- Optionally add extension fields later if needed.
- Find primitive roots for prime fields.
- Construct multiplicative subgroups of size `n` when `n | p - 1`.
- Construct cosets of such subgroups.

Required API:

```python
def is_prime(p: int) -> bool: ...
def factor(n: int) -> dict[int, int]: ...
def primitive_root(p: int) -> int: ...
def smooth_domain(p: int, n: int, coset_shift: int = 1) -> list[int]: ...
```

Acceptance:

- `smooth_domain(p, n)` rejects invalid `(p, n)` when `n` does not divide `p - 1`.
- Returned domain has exactly `n` distinct nonzero elements.
- Returned domain is closed under multiplication when `coset_shift = 1`.

### `rs_code.py`

Reed-Solomon encoding and code enumeration for tiny parameters.

Responsibilities:

- Evaluate degree `< k` polynomials on a domain.
- Enumerate all scalar RS codewords for tiny `(p, n, k)`.
- Fold scalar codewords into `GF(p)^m` symbols.
- Interleave `m` scalar codewords coordinatewise.
- Compute Hamming distance over scalar or tuple alphabets.

Required API:

```python
def eval_poly(coeffs: list[int], x: int, p: int) -> int: ...
def encode_rs(p: int, domain: list[int], coeffs: list[int]) -> tuple[int, ...]: ...
def enumerate_rs_codewords(p: int, domain: list[int], k: int) -> list[tuple[int, ...]]: ...
def fold_codeword(word: tuple[int, ...], m: int) -> tuple[tuple[int, ...], ...]: ...
def interleave_codewords(words: list[tuple[int, ...]]) -> tuple[tuple[int, ...], ...]: ...
def hamming_distance(a: tuple, b: tuple) -> int: ...
```

Acceptance:

- Number of scalar codewords equals `p^k`.
- Any two distinct scalar RS codewords have distance at least `n - k + 1`.
- Folding requires `m | n`.
- Interleaving requires all component codewords to have the same length.

### `exact_search.py`

Tiny-parameter exhaustive search.

Responsibilities:

- Enumerate all received words for very small alphabets and block lengths.
- Compute maximum list size at each integer radius.
- Record an example center attaining each maximum.
- Compare exact maxima against the volume lower bound and minimum-distance theorem.

Required API:

```python
def exact_max_list(codewords: list[tuple], alphabet: list, radius: int) -> tuple[int, tuple]: ...
def exact_profile(codewords: list[tuple], alphabet: list) -> list[dict]: ...
```

Acceptance:

- For scalar RS in fixed `n`, large-enough `p` toy regimes, exact profile matches the
  minimum-distance theorem:
  - below `(n - k + 1) / n`, list size is bounded by `binom(n, k)`;
  - at/above the minimum-weight boundary, the ball around zero contains many minimum-weight
    codewords.
- Search refuses cases whose received-word space is too large unless explicitly overridden.

### `structured_bad_centers.py`

Algebraic search for smooth-domain bad centers.

Responsibilities:

- Test received words of forms:
  - `r(x) = x^a`
  - `r(x) = x^a + lambda * x^b`
  - sparse sums of monomials
  - folded variants
- Count degree `< k` interpolants that agree with the center on many positions.
- Search for unusually large lists below the capacity candidate.

Required API:

```python
def monomial_center(p: int, domain: list[int], exponent: int) -> tuple[int, ...]: ...
def sparse_center(p: int, domain: list[int], terms: list[tuple[int, int]]) -> tuple[int, ...]: ...
def count_nearby_codewords(
    p: int,
    domain: list[int],
    k: int,
    center: tuple[int, ...],
    radius: int,
) -> int: ...
```

Acceptance:

- Produces reproducible search reports.
- Emits candidate counterexamples with full parameters and center definition.
- Distinguishes smooth-domain effects from random-domain effects.

## Experiments

### Experiment 1: Capacity Table

Generate capacity and volume-converse tables for:

```text
q_bits in {64, 128, 192, 256}
n in {2^10, 2^16, 2^20, 2^24}
m in {1, 2, 4, 8}
rho in {1/2, 1/4, 1/8, 1/16}
mode in {"folded", "interleaved"}
```

Output CSV columns:

```text
q_bits,n,m,rho,mode,N,Q_bits,entropy_target,delta_entropy,volume_grid_t,delta_volume_grid
```

### Experiment 2: Tiny Exact Enumeration

Use prime fields where `n | p - 1`, for example:

```text
p in {5, 7, 13, 17}
n in divisors of p - 1 with n <= 8
k in {1, ..., n - 1}
m in {1, 2} when valid
```

Compute exact max-list profiles and compare:

- Johnson radius.
- Minimum-distance boundary.
- Volume threshold.
- Capacity threshold.

### Experiment 3: Smooth vs Random Domains

For fixed `(p, n, k)`:

- Compare smooth subgroup domains.
- Compare cosets of smooth domains.
- Compare random subsets of size `n`.

Record whether structured domains show systematically larger lists near capacity.

### Experiment 4: Structured Bad Centers

Search exponents and sparse polynomial centers:

```text
a, b in [0, p - 2]
lambda in GF(p)^*
radius near floor(delta_cap * N)
```

Prioritize centers whose high-agreement subsets interpolate to degree `< k`.

## Proof Ledger

Maintain `docs/proof-ledger.md` with:

- Definitions and convention choices.
- Known upper bounds.
- Known lower bounds.
- Candidate lemmas.
- Failed proof attempts.
- Counterexample patterns.
- Exact theorem statements with required finite-size slack.

Initial candidate lemma:

```text
If H_{Q_m}(delta) <= 1 - R_m + log_{Q_m}(eps * q) / N_m - eta,
then every smooth-domain Reed-Solomon C_m satisfies
|Lambda(C_m, delta)| <= eps * q,
with explicit finite-size slack depending on eta, n, q, m, and rho.
```

This is the missing prize-level lemma. The volume bound already supplies the converse.

## Implementation Standards

- Keep experiments deterministic by default.
- Put all generated reports under `reports/`.
- Use log-space arithmetic for volume calculations.
- Refuse infeasible brute-force searches unless an override flag is supplied.
- Include parameters in every report filename.
- Keep folded and interleaved results separate in code, filenames, and tables.
- Add tests for mathematical invariants before scaling experiments.
- Use SageMath as an optional backend when available for extension fields, interpolation, and larger
  exact finite-field experiments. The core threshold and prime-field modules should remain pure
  Python so baseline reports can run without Sage.

## First Milestone

Create a minimal Python package with:

```text
src/rs_grand_list_decoding/rs_capacity_threshold.py
src/rs_grand_list_decoding/finite_field.py
src/rs_grand_list_decoding/rs_code.py
tests/test_capacity_threshold.py
tests/test_finite_field.py
tests/test_rs_code.py
```

Run:

```text
python -m pytest
python -m rs_grand_list_decoding.rs_capacity_threshold
```

The milestone is complete when:

- Capacity-table values match the pasted research.
- Tiny RS codewords satisfy the expected minimum distance.
- Smooth domains are constructed correctly for prime fields.
- Folded and interleaved modes are separately represented.

## Prize Strategy

1. Lock the convention by emailing the prize organizers about `C^{equiv m}` and the intended
   finite/asymptotic regime.
2. Build the numerical and exact-search harness.
3. Search for smooth-domain counterexamples below the capacity line.
4. If no counterexamples appear, focus proof work on the missing smooth-domain upper bound.
5. If counterexamples appear, generalize them into a lower-bound construction.
