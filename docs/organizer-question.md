# Organizer Convention Question

In the Grand List Decoding challenge, does "smooth domain" allow arbitrary multiplicative
subgroups/cosets whose size is `B`-smooth, e.g. `n = 900 = 2^2 * 3^2 * 5^2`, or is the intended
regime specifically power-of-two / 2-adic FFT domains?

This matters because for a subgroup `H` of size `n` and a subgroup `M <= H` of size `ell`, unions
of `M`-cosets give explicit lists of size `binom(n/ell, s/ell)` at radius `1 - s/n` whenever
`s > k` and `s - ell < k`. On power-of-two domains with dyadic rates this construction is tiny,
but on mixed-smooth domains it may exceed the `eps*q` budget.
