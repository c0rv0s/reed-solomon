# A Mixed-Smooth Lower-Bound Obstruction for the Proximity Prize Grand List-Decoding Challenge

This note records an explicit smooth-domain Reed-Solomon list-size obstruction under the natural
reading that a smooth evaluation domain may be a mixed-smooth multiplicative subgroup or coset.
If the intended convention is only power-of-two / 2-adic FFT domains, this example is instead a
request for convention clarification.

## Setup

Let `H <= F_p^*` be a multiplicative subgroup of size `n`, and let `RS[H,k]` be the scalar
Reed-Solomon code obtained by evaluating polynomials of degree `< k` on `H`.

Let `M <= H` be a subgroup of size `ell`. If `S` is a union of `r` cosets of `M`, write

```text
s = r * ell.
```

## Coset-Union Lemma

For a coset `aM`, the vanishing polynomial is

```text
prod_{alpha in aM} (x - alpha) = x^ell - a^ell.
```

Therefore, for a union `S` of `r` cosets,

```text
P_S(x) = prod_{alpha in S} (x - alpha)
       = prod_{i=1}^r (x^ell - y_i)
       = x^s + terms of degree at most s - ell.
```

For any received polynomial

```text
R(x) = x^s + h(x),    deg h < k,
```

define

```text
f_S(x) = R(x) - P_S(x).
```

If `s - ell < k`, then `deg f_S < k`. Also, for every `alpha in S`, `P_S(alpha)=0`, so
`f_S(alpha)=R(alpha)`. Thus `f_S` is a codeword agreeing with the received word on all `s` points
of `S`, i.e. within relative radius

```text
delta = 1 - s/n.
```

Distinct coset unions give distinct `P_S`, so the list size is at least

```text
binom(n/ell, s/ell).
```

The construction applies whenever:

```text
ell | n,  ell | s,  s > k,  s - ell < k.
```

## Concrete Witness

Use the prime

```text
p = 57896044618658832082471718862899876603311192244129634666496950740282297548801.
```

This prime has a short Proth certificate:

```text
p = 36028797018964425 * 2^200 + 1
36028797018964425 is odd
36028797018964425 < 2^200
7^((p-1)/2) = -1 mod p
```

By Proth's theorem, `p` is prime. The generated verification report also records:

```text
Proth certificate verified: true
p mod 900 = 1
log2(p) = 255.00000000000003
```

Since `p = 1 mod 900`, `F_p^*` contains a multiplicative subgroup `H` of size

```text
n = 900 = 2^2 * 3^2 * 5^2.
```

Take

```text
rho = 1/4
k = 225
ell = 4
M = n/ell = 225
r = 57
s = r*ell = 228
R(x) = x^228.
```

The conditions hold:

```text
s = 228 > 225 = k
s - ell = 224 < 225 = k.
```

Therefore, at radius

```text
delta = 1 - s/n = 1 - 228/900 = 0.7466666666666666,
```

the list size is at least

```text
binom(225,57)
  = 1230156107426602022802577569679062649087779656001965400
  > 2^179.68295852982317.
```

For `eps = 2^-128`, the allowed budget is

```text
eps * p = 2^-128 * p,
log2(eps * p) = 127.00000000000003.
```

Hence the list exceeds the budget by

```text
179.68295852982317 - 127.00000000000003 = 52.68295852982314 bits.
```

The entropy-capacity candidate computed by the repository for this concrete field size and
`n=900`, `rho=1/4`, `m=1`, folded/scalar mode is

```text
delta_entropy = 0.7473555627193669.
```

The obstruction radius satisfies

```text
0.7466666666666666 < 0.7473555627193669.
```

So the budget-violating list occurs below the entropy-capacity candidate.

## Conclusion

Under the public mixed-smooth-domain interpretation, this is a concrete counterexample to the
capacity-threshold answer for the grand list-decoding challenge: there is an explicit 256-bit prime
field, a smooth subgroup of size `900`, rate `1/4`, and a received word with list size exceeding
`eps * |F_p|` below the entropy-capacity candidate.

If the intended challenge only permits power-of-two / 2-adic FFT domains, then this example does
not disprove that narrower statement. It identifies the convention that must be made explicit.
