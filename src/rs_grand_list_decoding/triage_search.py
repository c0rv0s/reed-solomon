"""Sampled scale-and-triage search for structured RS bad centers."""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import math
import random
import re
from pathlib import Path
from typing import Any, Iterable

from rs_grand_list_decoding.exact_search import divisors
from rs_grand_list_decoding.finite_field import random_domain, smooth_domain
from rs_grand_list_decoding.structured_bad_centers import (
    add_score_metrics,
    lagrange_interpolate,
    monomial_center,
    poly_degree,
    sparse_center,
)
from rs_grand_list_decoding.rs_code import encode_rs, hamming_distance


def sampled_index_subsets(
    N: int,
    s: int,
    sample_budget: int,
    seed: int,
) -> tuple[list[tuple[int, ...]], bool]:
    """Return agreement subsets, exact if the full combination set fits the budget."""
    if not (0 <= s <= N):
        raise ValueError("s must be in [0, N]")
    if sample_budget <= 0:
        raise ValueError("sample_budget must be positive")

    total = math.comb(N, s)
    if total <= sample_budget:
        return list(itertools.combinations(range(N), s)), True

    rng = random.Random(seed)
    subsets: set[tuple[int, ...]] = set()
    max_attempts = sample_budget * 20
    attempts = 0
    while len(subsets) < sample_budget and attempts < max_attempts:
        subsets.add(tuple(sorted(rng.sample(range(N), s))))
        attempts += 1
    return sorted(subsets), False


def stable_seed(*parts: Any) -> int:
    """Return a deterministic 32-bit seed from structured values."""
    digest = hashlib.sha256("|".join(str(part) for part in parts).encode("utf-8")).digest()
    return int.from_bytes(digest[:4], "big")


def count_by_sampled_subsets(
    p: int,
    domain: list[int],
    k: int,
    center: tuple[int, ...],
    radius: int,
    sample_budget: int,
    seed: int,
) -> dict[str, Any]:
    """Find nearby codewords by sampled agreement-subset interpolation."""
    N = len(domain)
    if len(center) != N:
        raise ValueError("center length must match domain length")
    if radius < 0 or radius > N:
        raise ValueError("radius must be in [0, N]")
    s = N - radius
    if s < k:
        raise ValueError("sampled subset method requires agreement s >= k")

    subsets, exact = sampled_index_subsets(N, s, sample_budget=sample_budget, seed=seed)
    candidates: set[tuple[int, ...]] = set()
    for indices in subsets:
        xs = [domain[i] for i in indices]
        ys = [center[i] for i in indices]
        coeffs = lagrange_interpolate(p, xs, ys)
        if poly_degree(coeffs) >= k:
            continue
        codeword = encode_rs(p, domain, coeffs)
        if hamming_distance(center, codeword) <= radius:
            candidates.add(tuple(coeffs))
    return {
        "count_lower_bound": len(candidates),
        "subsets_checked": len(subsets),
        "subset_space": math.comb(N, s),
        "exact_subset_scan": exact,
    }


def _candidate_ks(n: int) -> list[int]:
    candidates = {n // 2, n // 4, n // 8}
    return sorted(k for k in candidates if 0 < k < n)


def _domain_specs(p: int, n: int, random_seeds: Iterable[int]) -> list[tuple[str, str, list[int]]]:
    specs = [("smooth", "subgroup", smooth_domain(p, n))]
    if n > 1:
        specs.append(("coset", "shift=2", smooth_domain(p, n, coset_shift=2)))
    for seed in random_seeds:
        specs.append(("random", f"seed={seed}", random_domain(p, n, seed=seed)))
    return specs


def _center_candidates(
    p: int,
    domain: list[int],
    exponent_limit: int,
    lambda_limit: int,
    constant_limit: int,
) -> Iterable[tuple[str, str, tuple[int, ...]]]:
    exponents = list(range(min(exponent_limit, p - 1)))
    lambdas = list(range(1, min(lambda_limit + 1, p)))
    constants = list(range(1, min(constant_limit + 1, p)))

    for a in exponents:
        yield "monomial", f"exponent={a}", monomial_center(p, domain, a)

    for a, b in itertools.combinations(exponents, 2):
        for lam in lambdas:
            yield "binomial", f"a={a};b={b};lambda={lam}", sparse_center(
                p, domain, [(1, a), (lam, b)]
            )
            for c in constants:
                yield "shifted_binomial", f"c={c};a={a};b={b};lambda={lam}", sparse_center(
                    p, domain, [(c, 0), (1, a), (lam, b)]
                )

    for a, b, c in itertools.combinations(exponents, 3):
        for lam1 in lambdas:
            for lam2 in lambdas:
                yield "trinomial", f"a={a};b={b};c={c};lambda1={lam1};lambda2={lam2}", sparse_center(
                    p, domain, [(1, a), (lam1, b), (lam2, c)]
                )


def _parse_parameters(parameters: str) -> dict[str, int]:
    return {key: int(value) for key, value in re.findall(r"([a-zA-Z0-9_]+)=(-?\d+)", parameters)}


def pattern_fields(center_type: str, parameters: str, n: int) -> dict[str, Any]:
    """Extract modular pattern fields from center parameters."""
    parsed = _parse_parameters(parameters)
    out: dict[str, Any] = {"pattern_center_type": center_type}
    for key in ("exponent", "a", "b", "c"):
        if key in parsed:
            out[f"{key}_mod_n"] = parsed[key] % n
    if "a" in parsed and "b" in parsed:
        out["b_minus_a_mod_n"] = (parsed["b"] - parsed["a"]) % n
    if "b" in parsed and "c" in parsed:
        out["c_minus_b_mod_n"] = (parsed["c"] - parsed["b"]) % n
    for key in ("lambda", "lambda1", "lambda2"):
        if key in parsed:
            out[key] = parsed[key]
    return out


def pattern_signature(row: dict[str, Any]) -> str:
    fields = pattern_fields(str(row["center_type"]), str(row["center_parameters"]), int(row["n"]))
    return ";".join(f"{key}={fields[key]}" for key in sorted(fields))


def _dedupe_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    seen: set[tuple] = set()
    deduped = []
    for row in rows:
        key = (
            row["p"],
            row["n"],
            row["k"],
            row["domain_type"],
            row["domain_label"],
            row["radius_label"],
            row["center_type"],
            tuple(row["center_signature"]),
        )
        if key in seen:
            continue
        seen.add(key)
        row = dict(row)
        row.pop("center_signature", None)
        deduped.append(row)
    return deduped


def generate_triage_rows(
    primes: Iterable[int] = (97, 193),
    max_n: int = 16,
    random_seed_count: int = 10,
    agreement_offsets: Iterable[int] = (1, 2, 3, 4),
    exponent_limit: int = 5,
    lambda_limit: int = 1,
    constant_limit: int = 1,
    sample_budget: int = 100,
) -> list[dict[str, Any]]:
    """Generate sampled scale-and-triage rows."""
    rows: list[dict[str, Any]] = []
    for p in primes:
        for n in divisors(p - 1):
            if n <= 1 or n > max_n:
                continue
            for k in _candidate_ks(n):
                for offset in agreement_offsets:
                    s = k + offset
                    if s > n:
                        continue
                    radius = n - s
                    radius_label = f"agreement-k+{offset}"
                    for domain_type, domain_label, domain in _domain_specs(
                        p, n, range(random_seed_count)
                    ):
                        for center_type, parameters, center in _center_candidates(
                            p,
                            domain,
                            exponent_limit=exponent_limit,
                            lambda_limit=lambda_limit,
                            constant_limit=constant_limit,
                        ):
                            result = count_by_sampled_subsets(
                                p,
                                domain,
                                k,
                                center,
                                radius,
                                sample_budget=sample_budget,
                                seed=stable_seed(p, n, k, radius_label, domain_label, parameters),
                            )
                            row = {
                                "p": p,
                                "n": n,
                                "k": k,
                                "rho": k / n,
                                "mode": "scalar",
                                "domain_type": domain_type,
                                "domain_label": domain_label,
                                "radius_label": radius_label,
                                "radius": radius,
                                "relative_radius": radius / n,
                                "center_type": center_type,
                                "center_parameters": parameters,
                                "count": result["count_lower_bound"],
                                "count_lower_bound": result["count_lower_bound"],
                                "subsets_checked": result["subsets_checked"],
                                "subset_space": result["subset_space"],
                                "exact_subset_scan": result["exact_subset_scan"],
                                "center_signature": center,
                                "method": "sampled_subset_interpolation",
                            }
                            add_score_metrics(row, p=p, N=n, k=k, radius=radius, m=1)
                            row.update(pattern_fields(center_type, parameters, n))
                            row["pattern_signature"] = pattern_signature(row)
                            rows.append(row)
    return sorted(
        _dedupe_rows(rows),
        key=lambda row: (
            bool(row["boundary_case"]),
            -float(row["generic_ratio"]),
            -int(row["count"]),
            row["p"],
            row["n"],
        ),
    )


def build_triage_comparison_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Compare smooth best rows with random best rows for each parameter family."""
    grouped: dict[tuple, list[dict[str, Any]]] = {}
    for row in rows:
        key = (
            row["p"],
            row["n"],
            row["k"],
            row["radius_label"],
            row["agreement_required"],
            row["center_type"],
        )
        grouped.setdefault(key, []).append(row)

    comparison_rows = []
    for key, group_rows in grouped.items():
        smooth_best = max(
            (int(row["count"]) for row in group_rows if row["domain_type"] == "smooth"),
            default=None,
        )
        random_best_by_domain: dict[str, int] = {}
        for row in group_rows:
            if row["domain_type"] == "random":
                label = str(row["domain_label"])
                random_best_by_domain[label] = max(random_best_by_domain.get(label, 0), int(row["count"]))
        if smooth_best is None or not random_best_by_domain:
            continue
        random_values = list(random_best_by_domain.values())
        random_mean = sum(random_values) / len(random_values)
        random_max = max(random_values)
        comparison_rows.append(
            {
                "p": key[0],
                "n": key[1],
                "k": key[2],
                "radius_label": key[3],
                "agreement_required": key[4],
                "center_type": key[5],
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
            -float(row["smooth_over_random_max"]),
            -float(row["smooth_over_random_mean"]),
            -int(row["smooth_best_count"]),
        ),
    )


def extract_pattern_rows(rows: list[dict[str, Any]], min_count: int = 2) -> list[dict[str, Any]]:
    """Group high rows by modular exponent pattern."""
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        if int(row["count"]) < min_count or bool(row["boundary_case"]):
            continue
        grouped.setdefault(str(row["pattern_signature"]), []).append(row)

    pattern_rows = []
    for signature, group_rows in grouped.items():
        best = max(group_rows, key=lambda row: (float(row["generic_ratio"]), int(row["count"])))
        pattern_rows.append(
            {
                "pattern_signature": signature,
                "occurrences": len(group_rows),
                "max_count": max(int(row["count"]) for row in group_rows),
                "max_generic_ratio": max(float(row["generic_ratio"]) for row in group_rows),
                "p_values": " ".join(str(p) for p in sorted({row["p"] for row in group_rows})),
                "n_values": " ".join(str(n) for n in sorted({row["n"] for row in group_rows})),
                "domain_types": " ".join(sorted({str(row["domain_type"]) for row in group_rows})),
                "best_parameters": best["center_parameters"],
            }
        )
    return sorted(
        pattern_rows,
        key=lambda row: (
            -int(row["occurrences"]),
            -float(row["max_generic_ratio"]),
            -int(row["max_count"]),
        ),
    )


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_triage_summary(path: Path, rows: list[dict[str, Any]], comparison_rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    top_rows = rows[:40]
    lines = [
        "# Triage Search Summary",
        "",
        "Sampled subset-interpolation sweep for larger prime fields. Counts are lower bounds unless `exact_subset_scan` is true.",
        "",
        f"- Rows: `{len(rows)}`",
        f"- Comparison rows: `{len(comparison_rows)}`",
        "",
        "| count | clipped ratio | raw ratio | p | n | k | s | domain | type | parameters | exact subsets |",
        "|---:|---:|---:|---:|---:|---:|---:|---|---|---|---|",
    ]
    for row in top_rows:
        lines.append(
            "| {count} | {generic_ratio:.3f} | {generic_ratio_raw:.3f} | {p} | {n} | {k} | {agreement_required} | {domain_type}:{domain_label} | {center_type} | `{center_parameters}` | {exact_subset_scan} |".format(
                **row
            )
        )
    path.write_text("\n".join(lines) + "\n")


def write_comparison_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Triage Domain Comparison Summary",
        "",
        "Smooth subgroup best sampled counts compared with random-domain best sampled counts.",
        "",
        f"- Rows: `{len(rows)}`",
        "",
        "| smooth/random max | smooth/random mean | smooth best | random max | random mean | p | n | k | s | type |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in rows[:40]:
        lines.append(
            "| {smooth_over_random_max:.3f} | {smooth_over_random_mean:.3f} | {smooth_best_count} | {random_max_best_count} | {random_mean_best_count:.3f} | {p} | {n} | {k} | {agreement_required} | {center_type} |".format(
                **row
            )
        )
    path.write_text("\n".join(lines) + "\n")


def write_patterns_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Triage Pattern Summary",
        "",
        "High-count non-boundary rows grouped by modular exponent pattern.",
        "",
        f"- Rows: `{len(rows)}`",
        "",
        "| occurrences | max ratio | max count | p values | n values | domains | pattern | example |",
        "|---:|---:|---:|---|---|---|---|---|",
    ]
    for row in rows[:40]:
        lines.append(
            "| {occurrences} | {max_generic_ratio:.3f} | {max_count} | {p_values} | {n_values} | {domain_types} | `{pattern_signature}` | `{best_parameters}` |".format(
                **row
            )
        )
    path.write_text("\n".join(lines) + "\n")


TRIAGE_FIELDS = [
    "p",
    "n",
    "k",
    "rho",
    "mode",
    "domain_type",
    "domain_label",
    "radius_label",
    "radius",
    "relative_radius",
    "agreement_required",
    "subset_count",
    "generic_expected",
    "generic_ratio",
    "generic_ratio_raw",
    "boundary_ratio",
    "boundary_case",
    "center_type",
    "center_parameters",
    "count",
    "count_lower_bound",
    "subsets_checked",
    "subset_space",
    "exact_subset_scan",
    "pattern_signature",
    "method",
]


COMPARISON_FIELDS = [
    "p",
    "n",
    "k",
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


PATTERN_FIELDS = [
    "pattern_signature",
    "occurrences",
    "max_count",
    "max_generic_ratio",
    "p_values",
    "n_values",
    "domain_types",
    "best_parameters",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--comparison-csv", type=Path, required=True)
    parser.add_argument("--comparison-summary", type=Path, required=True)
    parser.add_argument("--patterns-csv", type=Path, required=True)
    parser.add_argument("--patterns-summary", type=Path, required=True)
    parser.add_argument("--random-seeds", type=int, default=10)
    parser.add_argument("--sample-budget", type=int, default=100)
    parser.add_argument("--max-n", type=int, default=16)
    parser.add_argument("--exponent-limit", type=int, default=5)
    args = parser.parse_args()

    rows = generate_triage_rows(
        random_seed_count=args.random_seeds,
        sample_budget=args.sample_budget,
        max_n=args.max_n,
        exponent_limit=args.exponent_limit,
    )
    comparison_rows = build_triage_comparison_rows(rows)
    pattern_rows = extract_pattern_rows(rows)

    write_csv(args.csv, rows, TRIAGE_FIELDS)
    write_csv(args.comparison_csv, comparison_rows, COMPARISON_FIELDS)
    write_csv(args.patterns_csv, pattern_rows, PATTERN_FIELDS)
    write_triage_summary(args.summary, rows, comparison_rows)
    write_comparison_summary(args.comparison_summary, comparison_rows)
    write_patterns_summary(args.patterns_summary, pattern_rows)


if __name__ == "__main__":
    main()
