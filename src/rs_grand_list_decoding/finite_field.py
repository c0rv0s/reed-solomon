"""Small prime-field utilities for Reed-Solomon experiments."""

from __future__ import annotations

import random


def is_prime(p: int) -> bool:
    """Return whether p is prime, using trial division suitable for experiments."""
    if p < 2:
        return False
    if p in (2, 3):
        return True
    if p % 2 == 0 or p % 3 == 0:
        return False
    d = 5
    while d * d <= p:
        if p % d == 0 or p % (d + 2) == 0:
            return False
        d += 6
    return True


def factor(n: int) -> dict[int, int]:
    """Return the prime factorization of a positive integer."""
    if n <= 0:
        raise ValueError("n must be positive")
    factors: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d = 3 if d == 2 else d + 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def primitive_root(p: int) -> int:
    """Return a primitive root modulo a prime p."""
    if not is_prime(p):
        raise ValueError("p must be prime")
    if p == 2:
        return 1

    prime_factors = factor(p - 1).keys()
    for g in range(2, p):
        if all(pow(g, (p - 1) // r, p) != 1 for r in prime_factors):
            return g
    raise RuntimeError(f"no primitive root found for prime {p}")


def smooth_domain(p: int, n: int, coset_shift: int = 1) -> list[int]:
    """Return a multiplicative subgroup or coset of size n in GF(p)^*."""
    if not is_prime(p):
        raise ValueError("p must be prime")
    if n <= 0:
        raise ValueError("n must be positive")
    if (p - 1) % n != 0:
        raise ValueError("n must divide p - 1")

    shift = coset_shift % p
    if shift == 0:
        raise ValueError("coset_shift must be nonzero modulo p")

    g = primitive_root(p)
    subgroup_generator = pow(g, (p - 1) // n, p)
    domain = [(shift * pow(subgroup_generator, i, p)) % p for i in range(n)]
    if len(set(domain)) != n:
        raise RuntimeError("constructed domain has duplicate elements")
    return domain


def random_domain(p: int, n: int, seed: int) -> list[int]:
    """Return a deterministic random size-n subset of GF(p)^*."""
    if not is_prime(p):
        raise ValueError("p must be prime")
    if n <= 0:
        raise ValueError("n must be positive")
    if n > p - 1:
        raise ValueError("n cannot exceed p - 1")
    rng = random.Random(seed)
    return rng.sample(range(1, p), n)
