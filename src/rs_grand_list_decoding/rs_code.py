"""Reed-Solomon encoding utilities for tiny exact experiments."""

from __future__ import annotations

from itertools import product


def eval_poly(coeffs: list[int], x: int, p: int) -> int:
    """Evaluate sum_i coeffs[i] * x^i modulo p."""
    if p <= 1:
        raise ValueError("p must be at least 2")
    value = 0
    for coeff in reversed(coeffs):
        value = (value * x + coeff) % p
    return value


def encode_rs(p: int, domain: list[int], coeffs: list[int]) -> tuple[int, ...]:
    """Evaluate a degree-<k polynomial on the given domain."""
    return tuple(eval_poly(coeffs, x, p) for x in domain)


def enumerate_rs_codewords(p: int, domain: list[int], k: int) -> list[tuple[int, ...]]:
    """Enumerate all scalar RS codewords for tiny parameters."""
    if k < 0:
        raise ValueError("k must be nonnegative")
    return [encode_rs(p, domain, list(coeffs)) for coeffs in product(range(p), repeat=k)]


def fold_codeword(word: tuple[int, ...], m: int) -> tuple[tuple[int, ...], ...]:
    """Group a scalar word into consecutive m-tuples."""
    if m <= 0:
        raise ValueError("m must be positive")
    if len(word) % m != 0:
        raise ValueError("folding requires m to divide the word length")
    return tuple(tuple(word[i : i + m]) for i in range(0, len(word), m))


def interleave_codewords(words: list[tuple[int, ...]]) -> tuple[tuple[int, ...], ...]:
    """Group m scalar codewords coordinatewise into tuple alphabet symbols."""
    if not words:
        raise ValueError("at least one word is required")
    n = len(words[0])
    if any(len(word) != n for word in words):
        raise ValueError("all words must have the same length")
    return tuple(tuple(word[i] for word in words) for i in range(n))


def hamming_distance(a: tuple, b: tuple) -> int:
    """Return Hamming distance between equal-length words."""
    if len(a) != len(b):
        raise ValueError("words must have equal length")
    return sum(x != y for x, y in zip(a, b, strict=True))

