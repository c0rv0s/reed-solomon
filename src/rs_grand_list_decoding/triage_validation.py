"""Validation reruns for top sampled triage candidates."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Any, Iterable

from rs_grand_list_decoding.finite_field import random_domain, smooth_domain
from rs_grand_list_decoding.structured_bad_centers import monomial_center, sparse_center
from rs_grand_list_decoding.triage_search import (
    _parse_parameters,
    count_by_sampled_subsets,
    stable_seed,
)


def parse_int_list(text: str) -> list[int]:
    """Parse comma/space-separated positive integers."""
    values = [part for chunk in text.split(",") for part in chunk.split()]
    parsed = [int(value) for value in values if value]
    if not parsed or any(value <= 0 for value in parsed):
        raise ValueError("integer list must contain positive values")
    return sorted(set(parsed))


def _int(row: dict[str, Any], key: str) -> int:
    return int(row[key])


def _float(row: dict[str, Any], key: str, default: float = 0.0) -> float:
    value = row.get(key, "")
    if value == "":
        return default
    return float(value)


def _bool(row: dict[str, Any], key: str, default: bool = False) -> bool:
    value = row.get(key, "")
    if value == "":
        return default
    return str(value).lower() == "true"


def center_from_row(p: int, domain: list[int], row: dict[str, Any]) -> tuple[int, ...]:
    """Rebuild the structured center encoded by a triage CSV row."""
    center_type = str(row["center_type"])
    parsed = _parse_parameters(str(row["center_parameters"]))
    if center_type == "monomial":
        return monomial_center(p, domain, parsed["exponent"])
    if center_type == "binomial":
        return sparse_center(p, domain, [(1, parsed["a"]), (parsed["lambda"], parsed["b"])])
    if center_type == "shifted_binomial":
        return sparse_center(
            p,
            domain,
            [(parsed["c"], 0), (1, parsed["a"]), (parsed["lambda"], parsed["b"])],
        )
    if center_type == "trinomial":
        return sparse_center(
            p,
            domain,
            [(1, parsed["a"]), (parsed["lambda1"], parsed["b"]), (parsed["lambda2"], parsed["c"])],
        )
    raise ValueError(f"unsupported center_type={center_type!r}")


def load_top_candidates(
    path: Path,
    top_limit: int,
    min_count: int = 2,
    smooth_only: bool = True,
    non_boundary: bool = True,
) -> list[dict[str, Any]]:
    """Load and rank top triage rows from a CSV report."""
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))

    candidates = []
    seen: set[tuple[Any, ...]] = set()
    for row in rows:
        if smooth_only and row.get("domain_type") != "smooth":
            continue
        if non_boundary and _bool(row, "boundary_case"):
            continue
        if _int(row, "count") < min_count:
            continue
        key = (
            row["p"],
            row["n"],
            row["k"],
            row["radius"],
            row["center_type"],
            row["center_parameters"],
        )
        if key in seen:
            continue
        seen.add(key)
        candidates.append(row)

    def sort_key(row: dict[str, Any]) -> tuple[float, float, int, int, int, int]:
        sampled_ratio = _float(row, "sampled_generic_ratio", default=_float(row, "generic_ratio"))
        raw_ratio = _float(
            row,
            "full_generic_ratio_lower_bound",
            default=_float(row, "generic_ratio_raw", default=_float(row, "generic_ratio")),
        )
        return (
            -sampled_ratio,
            -raw_ratio,
            -_int(row, "count"),
            _int(row, "p"),
            _int(row, "n"),
            _int(row, "k"),
        )

    return sorted(candidates, key=sort_key)[:top_limit]


def _validated_sample_budget(subset_space: int, sample_budget: int, exact_threshold: int) -> int:
    if subset_space <= exact_threshold:
        return subset_space
    return sample_budget


def _count_candidate(
    p: int,
    domain: list[int],
    row: dict[str, Any],
    sample_budget: int,
    exact_threshold: int,
    seed: int,
) -> dict[str, Any]:
    n = _int(row, "n")
    k = _int(row, "k")
    radius = _int(row, "radius")
    agreement = n - radius
    subset_space = math.comb(n, agreement)
    effective_budget = _validated_sample_budget(subset_space, sample_budget, exact_threshold)
    center = center_from_row(p, domain, row)
    return count_by_sampled_subsets(
        p,
        domain,
        k,
        center,
        radius,
        sample_budget=effective_budget,
        seed=seed,
    )


def validate_candidate(
    row: dict[str, Any],
    budgets: Iterable[int],
    random_seed_count: int,
    exact_threshold: int,
) -> dict[str, Any]:
    """Validate one smooth-domain triage candidate against random domains."""
    budgets = sorted(set(int(budget) for budget in budgets))
    p = _int(row, "p")
    n = _int(row, "n")
    k = _int(row, "k")
    radius = _int(row, "radius")
    agreement = n - radius
    subset_space = math.comb(n, agreement)
    smooth_domain_values = smooth_domain(p, n)

    smooth_counts: dict[int, int] = {}
    smooth_subsets: dict[int, int] = {}
    smooth_exact: dict[int, bool] = {}
    for budget in budgets:
        result = _count_candidate(
            p,
            smooth_domain_values,
            row,
            sample_budget=budget,
            exact_threshold=exact_threshold,
            seed=stable_seed(
                "validation-smooth",
                p,
                n,
                k,
                radius,
                row["center_type"],
                row["center_parameters"],
                budget,
            ),
        )
        smooth_counts[budget] = int(result["count_lower_bound"])
        smooth_subsets[budget] = int(result["subsets_checked"])
        smooth_exact[budget] = bool(result["exact_subset_scan"])

    validation_budget = max(budgets)
    random_counts = []
    for random_seed in range(random_seed_count):
        domain = random_domain(p, n, seed=random_seed)
        result = _count_candidate(
            p,
            domain,
            row,
            sample_budget=validation_budget,
            exact_threshold=exact_threshold,
            seed=stable_seed(
                "validation-random",
                p,
                n,
                k,
                radius,
                row["center_type"],
                row["center_parameters"],
                random_seed,
                validation_budget,
            ),
        )
        random_counts.append(int(result["count_lower_bound"]))

    smooth_best = max(smooth_counts.values())
    random_mean = sum(random_counts) / len(random_counts) if random_counts else 0.0
    random_max = max(random_counts) if random_counts else 0
    p_value = (
        sum(1 for count in random_counts if count >= smooth_best) / len(random_counts)
        if random_counts
        else math.nan
    )
    effective_validation_budget = _validated_sample_budget(
        subset_space,
        validation_budget,
        exact_threshold,
    )
    out: dict[str, Any] = {
        "p": p,
        "n": n,
        "k": k,
        "agreement_required": agreement,
        "radius": radius,
        "subset_space": subset_space,
        "exact_subset_scan": subset_space <= exact_threshold,
        "validated_sample_budget": effective_validation_budget,
        "validated_sample_fraction": effective_validation_budget / subset_space,
        "center_type": row["center_type"],
        "center_parameters": row["center_parameters"],
        "input_count": _int(row, "count"),
        "input_sampled_generic_ratio": _float(
            row,
            "sampled_generic_ratio",
            default=_float(row, "generic_ratio"),
        ),
        "input_full_generic_ratio_lower_bound": _float(
            row,
            "full_generic_ratio_lower_bound",
            default=_float(row, "generic_ratio_raw", default=_float(row, "generic_ratio")),
        ),
        "low_k_fiber_predicted_count": row.get("low_k_fiber_predicted_count", ""),
        "low_k_fiber_explained": row.get("low_k_fiber_explained", ""),
        "smooth_best_count": smooth_best,
        "random_seed_count": random_seed_count,
        "random_mean_count": random_mean,
        "random_max_count": random_max,
        "empirical_p_value": p_value,
        "random_counts_sample": " ".join(str(count) for count in random_counts[:20]),
    }
    for budget in budgets:
        out[f"smooth_count_budget_{budget}"] = smooth_counts[budget]
        out[f"smooth_subsets_budget_{budget}"] = smooth_subsets[budget]
        out[f"smooth_exact_budget_{budget}"] = smooth_exact[budget]
    return out


def validate_candidates(
    candidates: Iterable[dict[str, Any]],
    budgets: Iterable[int],
    random_seed_count: int,
    exact_threshold: int,
) -> list[dict[str, Any]]:
    """Validate a sequence of candidate rows."""
    return [
        validate_candidate(
            row,
            budgets=budgets,
            random_seed_count=random_seed_count,
            exact_threshold=exact_threshold,
        )
        for row in candidates
    ]


def validation_fields(budgets: Iterable[int]) -> list[str]:
    """Return CSV field order for validation output."""
    fields = [
        "p",
        "n",
        "k",
        "agreement_required",
        "radius",
        "subset_space",
        "exact_subset_scan",
        "validated_sample_budget",
        "validated_sample_fraction",
        "center_type",
        "center_parameters",
        "input_count",
        "input_sampled_generic_ratio",
        "input_full_generic_ratio_lower_bound",
        "low_k_fiber_predicted_count",
        "low_k_fiber_explained",
    ]
    for budget in sorted(set(int(budget) for budget in budgets)):
        fields.extend(
            [
                f"smooth_count_budget_{budget}",
                f"smooth_subsets_budget_{budget}",
                f"smooth_exact_budget_{budget}",
            ]
        )
    fields.extend(
        [
            "smooth_best_count",
            "random_seed_count",
            "random_mean_count",
            "random_max_count",
            "empirical_p_value",
            "random_counts_sample",
        ]
    )
    return fields


def write_validation_csv(path: Path, rows: list[dict[str, Any]], budgets: Iterable[int]) -> None:
    """Write validation rows to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=validation_fields(budgets), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_validation_summary(
    path: Path,
    rows: list[dict[str, Any]],
    budgets: Iterable[int],
    input_path: Path,
) -> None:
    """Write compact Markdown summary for validation reruns."""
    path.parent.mkdir(parents=True, exist_ok=True)
    budget_text = ", ".join(str(budget) for budget in sorted(set(int(b) for b in budgets)))
    lines = [
        "# Triage Validation Summary",
        "",
        "Reruns top smooth-domain triage candidates at larger subset budgets, then compares the",
        "same center family against deterministic random domains with matching parameters.",
        "",
        f"- Input: `{input_path}`",
        f"- Rows: `{len(rows)}`",
        f"- Smooth budgets: `{budget_text}`",
        "",
        "| p-value | smooth best | random max | random mean | p | n | k | s | type | parameters | low-k explained |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|",
    ]
    for row in sorted(rows, key=lambda item: (float(item["empirical_p_value"]), -int(item["smooth_best_count"]))):
        lines.append(
            "| {empirical_p_value:.3f} | {smooth_best_count} | {random_max_count} | {random_mean_count:.3f} | {p} | {n} | {k} | {agreement_required} | {center_type} | `{center_parameters}` | {low_k_fiber_explained} |".format(
                **row
            )
        )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--budgets", default="100,1000,10000")
    parser.add_argument("--top-limit", type=int, default=10)
    parser.add_argument("--min-count", type=int, default=2)
    parser.add_argument("--exact-threshold", type=int, default=1_000_000)
    parser.add_argument("--random-seed-count", type=int, default=100)
    args = parser.parse_args()

    budgets = parse_int_list(args.budgets)
    candidates = load_top_candidates(
        args.input,
        top_limit=args.top_limit,
        min_count=args.min_count,
        smooth_only=True,
        non_boundary=True,
    )
    rows = validate_candidates(
        candidates,
        budgets=budgets,
        random_seed_count=args.random_seed_count,
        exact_threshold=args.exact_threshold,
    )
    write_validation_csv(args.csv, rows, budgets)
    write_validation_summary(args.summary, rows, budgets, args.input)


if __name__ == "__main__":
    main()
