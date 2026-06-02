"""Structured bad-center searches for smooth-domain Reed-Solomon codes."""

from __future__ import annotations

import argparse
import csv
import itertools
import math
from pathlib import Path
from typing import Any, Iterable

from rs_grand_list_decoding.exact_search import DEFAULT_MAX_CODEWORDS, divisors
from rs_grand_list_decoding.finite_field import random_domain, smooth_domain
from rs_grand_list_decoding.rs_code import (
    encode_rs,
    enumerate_rs_codewords,
    fold_codeword,
    hamming_distance,
)


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


def count_nearby_folded_codewords_bruteforce(
    p: int,
    domain: list[int],
    k: int,
    center: tuple[int, ...],
    radius: int,
    m: int,
    max_codewords: int = DEFAULT_MAX_CODEWORDS,
) -> int:
    """Count folded RS codewords within tuple-symbol radius by brute force."""
    if len(domain) % m != 0:
        raise ValueError("folded mode requires m | len(domain)")
    if p**k > max_codewords:
        raise ValueError(
            f"folded search would enumerate {p**k} codewords, "
            f"exceeding max_codewords={max_codewords}"
        )
    folded_center = fold_codeword(center, m)
    folded_codewords = [
        fold_codeword(codeword, m) for codeword in enumerate_rs_codewords(p, domain, k)
    ]
    return sum(hamming_distance(folded_center, codeword) <= radius for codeword in folded_codewords)


def baseline_metrics(p: int, N: int, k: int, radius: int, m: int = 1) -> dict[str, Any]:
    """Return generic baseline metrics for a symbol-agreement search row."""
    if radius < 0 or radius > N:
        raise ValueError("radius must be in [0, N]")
    agreement_required = N - radius
    subset_count = math.comb(N, agreement_required)
    constraints = m * agreement_required
    generic_expected = subset_count * (p ** (k - constraints))
    return {
        "agreement_required": agreement_required,
        "subset_count": subset_count,
        "generic_expected": generic_expected,
        "boundary_case": (agreement_required - 1) * m < k <= agreement_required * m,
    }


def add_score_metrics(row: dict[str, Any], p: int, N: int, k: int, radius: int, m: int) -> None:
    """Mutate a row with baseline-derived score metrics."""
    metrics = baseline_metrics(p, N, k, radius, m=m)
    row.update(metrics)
    row["generic_ratio"] = row["count"] / max(1.0, float(metrics["generic_expected"]))
    row["boundary_ratio"] = row["count"] / max(1, int(metrics["subset_count"]))


def _deduped_scored_rows(
    p: int,
    domain: list[int],
    k: int,
    radius: int,
    center_type: str,
    candidates: Iterable[tuple[str, tuple[int, ...]]],
    domain_type: str,
    domain_label: str,
    mode: str = "scalar",
    m: int = 1,
    max_folded_codewords: int = DEFAULT_MAX_CODEWORDS,
) -> list[dict[str, Any]]:
    """Score centers once per distinct restricted center signature."""
    n = len(domain)
    N = n // m if mode == "folded" else n
    groups: dict[tuple, dict[str, Any]] = {}

    for parameter, center in candidates:
        signature = fold_codeword(center, m) if mode == "folded" else center
        if signature in groups:
            group = groups[signature]
            group["duplicate_parameter_count"] += 1
            if len(group["all_parameters_sample"]) < 5:
                group["all_parameters_sample"].append(parameter)
            continue

        if mode == "folded":
            count = count_nearby_folded_codewords_bruteforce(
                p,
                domain,
                k,
                center,
                radius,
                m=m,
                max_codewords=max_folded_codewords,
            )
        elif mode == "scalar":
            count = count_nearby_codewords_by_subsets(p, domain, k, center, radius)
        else:
            raise ValueError("mode must be 'scalar' or 'folded'")

        row = {
            "p": p,
            "n": n,
            "k": k,
            "rho": k / n,
            "mode": mode,
            "m": m,
            "N": N,
            "domain_type": domain_type,
            "domain_label": domain_label,
            "radius": radius,
            "relative_radius": radius / N,
            "center_type": center_type,
            "center_parameters": parameter,
            "count": count,
            "duplicate_parameter_count": 1,
            "all_parameters_sample": [parameter],
            "method": "folded_bruteforce" if mode == "folded" else "subset_interpolation",
        }
        add_score_metrics(row, p=p, N=N, k=k, radius=radius, m=m)
        groups[signature] = row

    rows = []
    for row in groups.values():
        row["all_parameters_sample"] = " | ".join(row["all_parameters_sample"])
        rows.append(row)
    return sorted(
        rows,
        key=lambda row: (
            bool(row["boundary_case"]),
            -float(row["generic_ratio"]),
            -int(row["count"]),
            str(row["center_parameters"]),
        ),
    )


def search_monomial_centers(
    p: int,
    n: int,
    k: int,
    radius: int,
    exponent_range: Iterable[int] | None = None,
    domain: list[int] | None = None,
    domain_type: str = "smooth",
    domain_label: str = "subgroup",
    mode: str = "scalar",
    m: int = 1,
) -> list[dict[str, Any]]:
    """Search centers r(x)=x^a and return rows sorted by count descending."""
    domain = domain if domain is not None else smooth_domain(p, n)
    exponents = list(exponent_range if exponent_range is not None else range(p - 1))
    candidates = ((f"exponent={exponent}", monomial_center(p, domain, exponent)) for exponent in exponents)
    return _deduped_scored_rows(
        p,
        domain,
        k,
        radius,
        center_type="monomial",
        candidates=candidates,
        domain_type=domain_type,
        domain_label=domain_label,
        mode=mode,
        m=m,
    )


def search_binomial_centers(
    p: int,
    n: int,
    k: int,
    radius: int,
    exponent_range: Iterable[int] | None = None,
    lambdas: Iterable[int] | None = None,
    domain: list[int] | None = None,
    domain_type: str = "smooth",
    domain_label: str = "subgroup",
    mode: str = "scalar",
    m: int = 1,
) -> list[dict[str, Any]]:
    """Search centers r(x)=x^a+lambda*x^b and return rows sorted by count."""
    domain = domain if domain is not None else smooth_domain(p, n)
    exponents = list(exponent_range if exponent_range is not None else range(p - 1))
    lambda_values = list(lambdas if lambdas is not None else range(1, p))
    candidates = (
        (
            f"a={a};b={b};lambda={lam}",
            sparse_center(p, domain, [(1, a), (lam, b)]),
        )
        for a, b in itertools.combinations(exponents, 2)
        for lam in lambda_values
    )
    return _deduped_scored_rows(
        p,
        domain,
        k,
        radius,
        center_type="binomial",
        candidates=candidates,
        domain_type=domain_type,
        domain_label=domain_label,
        mode=mode,
        m=m,
    )


def _candidate_ks(n: int) -> list[int]:
    candidates = {n // 2, n // 4}
    return sorted(k for k in candidates if 0 < k < n)


def _agreement_radii(n: int, k: int, mode: str, m: int, offsets: Iterable[int]) -> dict[str, int]:
    N = n // m if mode == "folded" else n
    base_agreement = math.ceil(k / m) if mode == "folded" else k
    radii = {}
    for offset in offsets:
        agreement = base_agreement + offset
        if agreement <= N:
            label = "agreement-k" if offset == 0 else f"agreement-k+{offset}"
            radii[label] = N - agreement
    return radii


def _domain_specs(
    p: int,
    n: int,
    random_seeds: Iterable[int],
    coset_shifts: Iterable[int],
) -> list[tuple[str, str, list[int]]]:
    seen: set[tuple[int, ...]] = set()
    specs: list[tuple[str, str, list[int]]] = []

    def add(domain_type: str, label: str, domain: list[int]) -> None:
        signature = tuple(sorted(domain))
        if signature in seen:
            return
        seen.add(signature)
        specs.append((domain_type, label, domain))

    add("smooth", "subgroup", smooth_domain(p, n))
    for shift in coset_shifts:
        if shift % p != 0:
            add("coset", f"shift={shift}", smooth_domain(p, n, coset_shift=shift))
    for seed in random_seeds:
        add("random", f"seed={seed}", random_domain(p, n, seed=seed))
    return specs


def _search_domain_centers(
    p: int,
    n: int,
    k: int,
    radius_label: str,
    radius: int,
    domain_type: str,
    domain_label: str,
    domain: list[int],
    exponent_limit: int,
    lambda_limit: int,
    top_k: int,
    mode: str,
    m: int,
) -> list[dict[str, Any]]:
    exponents = range(min(exponent_limit, p - 1))
    lambdas = range(1, min(lambda_limit + 1, p))
    monomial_rows = search_monomial_centers(
        p,
        n,
        k,
        radius,
        domain=domain,
        domain_type=domain_type,
        domain_label=domain_label,
        mode=mode,
        m=m,
    )
    binomial_rows = search_binomial_centers(
        p,
        n,
        k,
        radius,
        exponent_range=exponents,
        lambdas=lambdas,
        domain=domain,
        domain_type=domain_type,
        domain_label=domain_label,
        mode=mode,
        m=m,
    )
    rows = []
    for row in monomial_rows[:top_k] + binomial_rows[:top_k]:
        rows.append({**row, "radius_label": radius_label})
    return rows


def generate_structured_sweep_rows(
    primes: Iterable[int] = (17, 29, 41),
    max_n: int = 8,
    exponent_limit: int = 5,
    lambda_limit: int = 3,
    top_k: int = 3,
    random_seeds: Iterable[int] = (0, 1),
    coset_shifts: Iterable[int] = (2,),
    agreement_offsets: Iterable[int] = (0, 1, 2, 3),
    folded_m: int = 2,
    max_folded_codewords: int = 5_000,
) -> list[dict[str, Any]]:
    """Generate a structured-center sweep over smooth, coset, and random domains."""
    rows: list[dict[str, Any]] = []
    for p in primes:
        for n in divisors(p - 1):
            if n <= 1 or n > max_n:
                continue
            domain_specs = _domain_specs(p, n, random_seeds, coset_shifts)
            for k in _candidate_ks(n):
                for mode, m in (("scalar", 1), ("folded", folded_m)):
                    if mode == "folded":
                        if n % m != 0 or p**k > max_folded_codewords:
                            continue
                    for radius_label, radius in _agreement_radii(
                        n, k, mode=mode, m=m, offsets=agreement_offsets
                    ).items():
                        for domain_type, domain_label, domain in domain_specs:
                            rows.extend(
                                _search_domain_centers(
                                    p,
                                    n,
                                    k,
                                    radius_label,
                                    radius,
                                    domain_type,
                                    domain_label,
                                    domain,
                                    exponent_limit,
                                    lambda_limit,
                                    top_k,
                                    mode=mode,
                                    m=m,
                                )
                            )
    return rows


def build_domain_comparison_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Compare smooth-domain best counts against random-domain best counts."""
    grouped: dict[tuple, list[dict[str, Any]]] = {}
    for row in rows:
        key = (
            row["p"],
            row["n"],
            row["k"],
            row["mode"],
            row["m"],
            row["radius_label"],
            row["agreement_required"],
            row["center_type"],
        )
        grouped.setdefault(key, []).append(row)

    comparison_rows = []
    for key, group_rows in grouped.items():
        smooth_counts = [
            int(row["count"]) for row in group_rows if row["domain_type"] == "smooth"
        ]
        random_best_by_domain: dict[str, int] = {}
        for row in group_rows:
            if row["domain_type"] == "random":
                label = str(row["domain_label"])
                random_best_by_domain[label] = max(
                    random_best_by_domain.get(label, 0), int(row["count"])
                )
        if not smooth_counts or not random_best_by_domain:
            continue

        smooth_best = max(smooth_counts)
        random_values = list(random_best_by_domain.values())
        random_mean = sum(random_values) / len(random_values)
        random_max = max(random_values)
        comparison_rows.append(
            {
                "p": key[0],
                "n": key[1],
                "k": key[2],
                "mode": key[3],
                "m": key[4],
                "radius_label": key[5],
                "agreement_required": key[6],
                "center_type": key[7],
                "smooth_best_count": smooth_best,
                "random_mean_best_count": random_mean,
                "random_max_best_count": random_max,
                "smooth_over_random_mean": smooth_best / max(1.0, random_mean),
                "smooth_over_random_max": smooth_best / max(1, random_max),
                "random_domains": len(random_values),
            }
        )
    return sorted(
        comparison_rows,
        key=lambda row: (
            -float(row["smooth_over_random_mean"]),
            -int(row["smooth_best_count"]),
        ),
    )


def write_structured_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write structured-center search rows to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "p",
        "n",
        "k",
        "rho",
        "mode",
        "m",
        "N",
        "domain_type",
        "domain_label",
        "radius_label",
        "radius",
        "relative_radius",
        "agreement_required",
        "subset_count",
        "generic_expected",
        "generic_ratio",
        "boundary_ratio",
        "boundary_case",
        "center_type",
        "center_parameters",
        "count",
        "duplicate_parameter_count",
        "all_parameters_sample",
        "method",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_structured_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write a compact Markdown summary of the structured-center sweep."""
    path.parent.mkdir(parents=True, exist_ok=True)
    top_rows = sorted(
        rows,
        key=lambda row: (
            bool(row["boundary_case"]),
            -float(row["generic_ratio"]),
            -int(row["count"]),
            row["p"],
            row["n"],
        ),
    )[:50]
    lines = [
        "# Structured Center Summary",
        "",
        "Structured monomial/binomial center sweep with deduplication and baseline scoring.",
        "",
        f"- Rows: `{len(rows)}`",
        "- Default primes: `{17, 29, 41}`",
        "- Default smooth-domain limit: `n <= 8`",
        "- Domains: smooth subgroup, cosets, and deterministic random nonzero domains.",
        "",
        "| count | generic ratio | boundary | p | n | k | s | mode | domain | type | parameters | dupes |",
        "|---:|---:|---|---:|---:|---:|---:|---|---|---|---|---:|",
    ]
    for row in top_rows:
        lines.append(
            "| {count} | {generic_ratio:.3f} | {boundary_case} | {p} | {n} | {k} | {agreement_required} | {mode} | {domain_type}:{domain_label} | {center_type} | `{center_parameters}` | {duplicate_parameter_count} |".format(
                **row
            )
        )
    path.write_text("\n".join(lines) + "\n")


def write_domain_comparison_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write smooth-vs-random domain comparison rows to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "p",
        "n",
        "k",
        "mode",
        "m",
        "radius_label",
        "agreement_required",
        "center_type",
        "smooth_best_count",
        "random_mean_best_count",
        "random_max_best_count",
        "smooth_over_random_mean",
        "smooth_over_random_max",
        "random_domains",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_domain_comparison_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write a Markdown summary of smooth-vs-random comparison rows."""
    path.parent.mkdir(parents=True, exist_ok=True)
    top_rows = rows[:40]
    lines = [
        "# Domain Comparison Summary",
        "",
        "Smooth subgroup best counts compared with deterministic random-domain best counts.",
        "",
        f"- Rows: `{len(rows)}`",
        "",
        "| smooth/random mean | smooth/random max | smooth best | random mean | random max | p | n | k | s | mode | type |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for row in top_rows:
        lines.append(
            "| {smooth_over_random_mean:.3f} | {smooth_over_random_max:.3f} | {smooth_best_count} | {random_mean_best_count:.3f} | {random_max_best_count} | {p} | {n} | {k} | {agreement_required} | {mode} | {center_type} |".format(
                **row
            )
        )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--comparison-csv", type=Path)
    parser.add_argument("--comparison-summary", type=Path)
    parser.add_argument("--max-n", type=int, default=8)
    parser.add_argument("--exponent-limit", type=int, default=5)
    parser.add_argument("--lambda-limit", type=int, default=3)
    parser.add_argument("--top-k", type=int, default=3)
    parser.add_argument("--random-seeds", type=int, default=2)
    parser.add_argument("--max-folded-codewords", type=int, default=5_000)
    args = parser.parse_args()

    rows = generate_structured_sweep_rows(
        max_n=args.max_n,
        exponent_limit=args.exponent_limit,
        lambda_limit=args.lambda_limit,
        top_k=args.top_k,
        random_seeds=range(args.random_seeds),
        max_folded_codewords=args.max_folded_codewords,
    )
    write_structured_csv(args.csv, rows)
    write_structured_summary(args.summary, rows)
    if args.comparison_csv or args.comparison_summary:
        comparison_rows = build_domain_comparison_rows(rows)
        if args.comparison_csv:
            write_domain_comparison_csv(args.comparison_csv, comparison_rows)
        if args.comparison_summary:
            write_domain_comparison_summary(args.comparison_summary, comparison_rows)


if __name__ == "__main__":
    main()
