"""Structured bad-center searches for smooth-domain Reed-Solomon codes."""

from __future__ import annotations

import argparse
import csv
import itertools
import math
from pathlib import Path
from typing import Any, Iterable

from rs_grand_list_decoding.exact_search import divisors
from rs_grand_list_decoding.finite_field import smooth_domain
from rs_grand_list_decoding.rs_capacity_threshold import threshold_params
from rs_grand_list_decoding.rs_code import encode_rs, enumerate_rs_codewords, hamming_distance


def _trim(coeffs: list[int]) -> list[int]:
    while coeffs and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs


def monomial_center(p: int, domain: list[int], exponent: int) -> tuple[int, ...]:
    """Evaluate r(x) = x^exponent on the domain."""
    return tuple(pow(x, exponent, p) for x in domain)


def sparse_center(p: int, domain: list[int], terms: list[tuple[int, int]]) -> tuple[int, ...]:
    """Evaluate a sparse polynomial center with terms `(coefficient, exponent)`."""
    values = []
    for x in domain:
        value = 0
        for coeff, exponent in terms:
            value = (value + coeff * pow(x, exponent, p)) % p
        values.append(value)
    return tuple(values)


def poly_add_mod(p: int, a: list[int], b: list[int]) -> list[int]:
    """Add two low-to-high coefficient polynomials modulo p."""
    length = max(len(a), len(b))
    out = [0] * length
    for i in range(length):
        out[i] = ((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p
    return _trim(out)


def poly_mul_mod(p: int, a: list[int], b: list[int]) -> list[int]:
    """Multiply two low-to-high coefficient polynomials modulo p."""
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] = (out[i + j] + ai * bj) % p
    return _trim(out)


def poly_degree(coeffs: list[int]) -> int:
    """Return polynomial degree, or -1 for the zero polynomial."""
    for i in range(len(coeffs) - 1, -1, -1):
        if coeffs[i] != 0:
            return i
    return -1


def lagrange_interpolate(p: int, xs: list[int], ys: list[int]) -> list[int]:
    """Return low-to-high coefficients of the interpolant over GF(p)."""
    if len(xs) != len(ys):
        raise ValueError("xs and ys must have the same length")
    if len(set(x % p for x in xs)) != len(xs):
        raise ValueError("interpolation points must be distinct modulo p")

    result: list[int] = []
    for j, (xj, yj) in enumerate(zip(xs, ys, strict=True)):
        basis = [1]
        denom = 1
        for i, xi in enumerate(xs):
            if i == j:
                continue
            basis = poly_mul_mod(p, basis, [(-xi) % p, 1])
            denom = (denom * (xj - xi)) % p
        scale = (yj * pow(denom, -1, p)) % p
        result = poly_add_mod(p, result, [(scale * coeff) % p for coeff in basis])
    return _trim(result)


def count_nearby_codewords_bruteforce(
    p: int,
    domain: list[int],
    k: int,
    center: tuple[int, ...],
    radius: int,
) -> int:
    """Count degree-<k RS codewords within radius of a center by brute force."""
    if len(center) != len(domain):
        raise ValueError("center length must match domain length")
    if radius < 0 or radius > len(domain):
        raise ValueError("radius must be in [0, len(domain)]")
    codewords = enumerate_rs_codewords(p, domain, k)
    return sum(hamming_distance(center, codeword) <= radius for codeword in codewords)


def count_nearby_codewords_by_subsets(
    p: int,
    domain: list[int],
    k: int,
    center: tuple[int, ...],
    radius: int,
) -> int:
    """Count nearby codewords by interpolating agreement subsets and deduplicating."""
    N = len(domain)
    if len(center) != N:
        raise ValueError("center length must match domain length")
    if radius < 0 or radius > N:
        raise ValueError("radius must be in [0, len(domain)]")

    min_agreement = N - radius
    if min_agreement < k:
        return count_nearby_codewords_bruteforce(p, domain, k, center, radius)

    candidates: set[tuple[int, ...]] = set()
    for indices in itertools.combinations(range(N), min_agreement):
        xs = [domain[i] for i in indices]
        ys = [center[i] for i in indices]
        coeffs = lagrange_interpolate(p, xs, ys)
        if poly_degree(coeffs) >= k:
            continue
        codeword = encode_rs(p, domain, coeffs)
        if hamming_distance(center, codeword) <= radius:
            candidates.add(tuple(_trim(coeffs[:])))
    return len(candidates)


def search_monomial_centers(
    p: int,
    n: int,
    k: int,
    radius: int,
    exponent_range: Iterable[int] | None = None,
) -> list[dict[str, Any]]:
    """Search centers r(x)=x^a and return rows sorted by count descending."""
    domain = smooth_domain(p, n)
    exponents = list(exponent_range if exponent_range is not None else range(p - 1))
    rows = []
    for exponent in exponents:
        center = monomial_center(p, domain, exponent)
        count = count_nearby_codewords_by_subsets(p, domain, k, center, radius)
        rows.append(
            {
                "p": p,
                "n": n,
                "k": k,
                "rho": k / n,
                "radius": radius,
                "relative_radius": radius / n,
                "center_type": "monomial",
                "center_parameters": f"exponent={exponent}",
                "count": count,
                "method": "subset_interpolation",
            }
        )
    return sorted(rows, key=lambda row: (-row["count"], row["center_parameters"]))


def search_binomial_centers(
    p: int,
    n: int,
    k: int,
    radius: int,
    exponent_range: Iterable[int] | None = None,
    lambdas: Iterable[int] | None = None,
) -> list[dict[str, Any]]:
    """Search centers r(x)=x^a+lambda*x^b and return rows sorted by count."""
    domain = smooth_domain(p, n)
    exponents = list(exponent_range if exponent_range is not None else range(p - 1))
    lambda_values = list(lambdas if lambdas is not None else range(1, p))
    rows = []
    for a, b in itertools.combinations(exponents, 2):
        for lam in lambda_values:
            center = sparse_center(p, domain, [(1, a), (lam, b)])
            count = count_nearby_codewords_by_subsets(p, domain, k, center, radius)
            rows.append(
                {
                    "p": p,
                    "n": n,
                    "k": k,
                    "rho": k / n,
                    "radius": radius,
                    "relative_radius": radius / n,
                    "center_type": "binomial",
                    "center_parameters": f"a={a};b={b};lambda={lam}",
                    "count": count,
                    "method": "subset_interpolation",
                }
            )
    return sorted(rows, key=lambda row: (-row["count"], row["center_parameters"]))


def _candidate_ks(n: int) -> list[int]:
    candidates = {n // 2, n // 4}
    return sorted(k for k in candidates if 0 < k < n)


def _candidate_radii(p: int, n: int, k: int) -> dict[str, int]:
    rho = k / n
    d = n - k + 1
    cap = threshold_params(q_bits=math.log2(p), n=n, rho=rho, m=1, eps_bits=0.0)
    candidates = {
        "johnson": math.floor((1.0 - math.sqrt(rho)) * n),
        "capacity": math.floor(float(cap["delta_entropy"]) * n),
        "d-1": d - 1,
    }
    clipped = {label: max(0, min(n, radius)) for label, radius in candidates.items()}
    return {
        label: radius
        for label, radius in clipped.items()
        if n - radius >= k
    }


def generate_structured_sweep_rows(
    primes: Iterable[int] = (17, 29, 41),
    max_n: int = 12,
    exponent_limit: int = 6,
    lambda_limit: int = 4,
    top_k: int = 5,
) -> list[dict[str, Any]]:
    """Generate a first structured-center sweep over smooth subgroup domains."""
    rows: list[dict[str, Any]] = []
    for p in primes:
        for n in divisors(p - 1):
            if n <= 1 or n > max_n:
                continue
            for k in _candidate_ks(n):
                for radius_label, radius in _candidate_radii(p, n, k).items():
                    monomial_rows = search_monomial_centers(p, n, k, radius)
                    for row in monomial_rows[:top_k]:
                        rows.append({**row, "radius_label": radius_label})

                    exponents = range(min(exponent_limit, p - 1))
                    lambdas = range(1, min(lambda_limit + 1, p))
                    binomial_rows = search_binomial_centers(
                        p, n, k, radius, exponent_range=exponents, lambdas=lambdas
                    )
                    for row in binomial_rows[:top_k]:
                        rows.append({**row, "radius_label": radius_label})
    return rows


def write_structured_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write structured-center search rows to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "p",
        "n",
        "k",
        "rho",
        "radius_label",
        "radius",
        "relative_radius",
        "center_type",
        "center_parameters",
        "count",
        "method",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_structured_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write a compact Markdown summary of the structured-center sweep."""
    path.parent.mkdir(parents=True, exist_ok=True)
    top_rows = sorted(rows, key=lambda row: (-int(row["count"]), row["p"], row["n"]))[:40]
    lines = [
        "# Structured Center Summary",
        "",
        "First sweep over monomial and limited binomial centers on smooth subgroup domains.",
        "",
        f"- Rows: `{len(rows)}`",
        "- Default primes: `{17, 29, 41}`",
        "- Default smooth-domain limit: `n <= 12`",
        "",
        "| count | p | n | k | radius label | radius | type | parameters |",
        "|---:|---:|---:|---:|---|---:|---|---|",
    ]
    for row in top_rows:
        lines.append(
            "| {count} | {p} | {n} | {k} | {radius_label} | {radius} | {center_type} | `{center_parameters}` |".format(
                **row
            )
        )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--max-n", type=int, default=12)
    parser.add_argument("--exponent-limit", type=int, default=6)
    parser.add_argument("--lambda-limit", type=int, default=4)
    parser.add_argument("--top-k", type=int, default=5)
    args = parser.parse_args()

    rows = generate_structured_sweep_rows(
        max_n=args.max_n,
        exponent_limit=args.exponent_limit,
        lambda_limit=args.lambda_limit,
        top_k=args.top_k,
    )
    write_structured_csv(args.csv, rows)
    write_structured_summary(args.summary, rows)


if __name__ == "__main__":
    main()
