"""Explicit smooth-domain lower-bound constructions."""

from __future__ import annotations

import argparse
import csv
import itertools
import math
from pathlib import Path
from typing import Any, Iterable

from rs_grand_list_decoding.exact_search import divisors
from rs_grand_list_decoding.finite_field import primitive_root, smooth_domain
from rs_grand_list_decoding.rs_capacity_threshold import RATES, ln_comb, threshold_params
from rs_grand_list_decoding.rs_code import encode_rs
from rs_grand_list_decoding.structured_bad_centers import (
    count_nearby_codewords_by_subsets,
    poly_add_mod,
    poly_degree,
    poly_mul_mod,
)


def subgroup_generator_for_size(p: int, n: int) -> int:
    """Return a generator of the size-n subgroup of GF(p)^*."""
    if n <= 0:
        raise ValueError("n must be positive")
    if (p - 1) % n != 0:
        raise ValueError("n must divide p - 1")
    return pow(primitive_root(p), (p - 1) // n, p)


def subgroup_elements(p: int, n: int) -> list[int]:
    """Return the multiplicative subgroup H <= GF(p)^* of size n."""
    return smooth_domain(p, n)


def subgroup_of_subgroup(p: int, n: int, ell: int) -> list[int]:
    """Return the size-ell subgroup M <= H, where |H| = n."""
    if ell <= 0:
        raise ValueError("ell must be positive")
    if n % ell != 0:
        raise ValueError("ell must divide n")
    generator = subgroup_generator_for_size(p, n)
    sub_generator = pow(generator, n // ell, p)
    return [pow(sub_generator, i, p) for i in range(ell)]


def cosets_of_subgroup(p: int, H: list[int], M: list[int]) -> list[list[int]]:
    """Return all distinct cosets of M inside H."""
    h_set = set(H)
    m_set = set(M)
    if not m_set or not m_set.issubset(h_set):
        raise ValueError("M must be a nonempty subgroup of H")

    seen: set[int] = set()
    cosets = []
    for representative in H:
        if representative in seen:
            continue
        coset = sorted((representative * element) % p for element in M)
        if not set(coset).issubset(h_set):
            raise ValueError("coset left H")
        seen.update(coset)
        cosets.append(coset)
    return cosets


def _trim(coeffs: list[int]) -> list[int]:
    while coeffs and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs


def _sub_mod(p: int, a: list[int], b: list[int]) -> list[int]:
    return poly_add_mod(p, a, [(-coeff) % p for coeff in b])


def vanishing_polynomial(p: int, S: Iterable[int]) -> list[int]:
    """Return P_S(x) = prod_{a in S}(x-a), with low-to-high coefficients."""
    coeffs = [1]
    for element in S:
        coeffs = poly_mul_mod(p, coeffs, [(-element) % p, 1])
    return coeffs


def _lower_bound_conditions_hold(n: int, k: int, s: int, ell: int) -> bool:
    if n <= 0 or k <= 0 or s <= 0 or ell <= 0:
        return False
    if n % ell != 0 or s % ell != 0:
        return False
    if s > n or s <= k or s - ell >= k:
        return False
    return True


def coset_union_lower_bound_count(n: int, k: int, s: int, ell: int) -> int:
    """Return the coset-union lower-bound list size when the theorem applies."""
    if not _lower_bound_conditions_hold(n, k, s, ell):
        return 0
    return math.comb(n // ell, s // ell)


def _log2_coset_union_lower_bound(n: int, k: int, s: int, ell: int) -> float:
    if not _lower_bound_conditions_hold(n, k, s, ell):
        return -math.inf
    return ln_comb(n // ell, s // ell) / math.log(2.0)


def _center_coeffs(p: int, k: int, s: int, h_coeffs: list[int] | None) -> list[int]:
    h = [coeff % p for coeff in (h_coeffs or [])]
    if poly_degree(h[:]) >= k:
        raise ValueError("h_coeffs must have degree < k")
    coeffs = h[:] + [0] * max(0, s + 1 - len(h))
    coeffs[s] = (coeffs[s] + 1) % p
    return _trim(coeffs)


def construct_coset_union_codewords(
    p: int,
    n: int,
    k: int,
    s: int,
    ell: int,
    h_coeffs: list[int] | None = None,
    max_emit: int = 1000,
) -> list[dict[str, Any]]:
    """Emit codewords from the coset-union lower-bound construction."""
    if max_emit <= 0:
        raise ValueError("max_emit must be positive")
    lower_bound = coset_union_lower_bound_count(n, k, s, ell)
    if lower_bound == 0:
        raise ValueError("invalid coset-union parameters")

    H = subgroup_elements(p, n)
    M = subgroup_of_subgroup(p, n, ell)
    cosets = cosets_of_subgroup(p, H, M)
    cosets_to_choose = s // ell
    center_coeffs = _center_coeffs(p, k, s, h_coeffs)
    center_word = encode_rs(p, H, center_coeffs)
    positions = {element: index for index, element in enumerate(H)}
    emitted = []

    for selected_indices in itertools.islice(
        itertools.combinations(range(len(cosets)), cosets_to_choose),
        max_emit,
    ):
        agreement_set = sorted(
            element for index in selected_indices for element in cosets[index]
        )
        vanishing = vanishing_polynomial(p, agreement_set)
        polynomial = _sub_mod(p, center_coeffs, vanishing)
        degree = poly_degree(polynomial[:])
        if degree >= k:
            raise AssertionError("constructed polynomial has degree >= k")
        codeword = encode_rs(p, H, polynomial)
        selected_agreements = sum(
            codeword[positions[element]] == center_word[positions[element]]
            for element in agreement_set
        )
        if selected_agreements != s:
            raise AssertionError("constructed polynomial does not agree on S")
        total_agreement_count = sum(
            a == b for a, b in zip(codeword, center_word, strict=True)
        )
        emitted.append(
            {
                "selected_coset_indices": selected_indices,
                "selected_cosets": [cosets[index] for index in selected_indices],
                "agreement_set": agreement_set,
                "polynomial": polynomial,
                "degree": degree,
                "selected_agreement_count": selected_agreements,
                "total_agreement_count": total_agreement_count,
            }
        )
    return emitted


def verify_against_counting(
    p: int,
    n: int,
    k: int,
    s: int,
    ell: int,
    h_coeffs: list[int] | None = None,
) -> dict[str, Any]:
    """Compare the construction with exact subset counting for a small case."""
    lower_bound = coset_union_lower_bound_count(n, k, s, ell)
    constructed = construct_coset_union_codewords(
        p,
        n,
        k,
        s,
        ell,
        h_coeffs=h_coeffs,
        max_emit=max(1, lower_bound),
    )
    H = subgroup_elements(p, n)
    center_coeffs = _center_coeffs(p, k, s, h_coeffs)
    exact_count = count_nearby_codewords_by_subsets(
        p,
        H,
        k,
        tuple(encode_rs(p, H, center_coeffs)),
        radius=n - s,
    )
    constructed_polynomials = {tuple(row["polynomial"]) for row in constructed}
    return {
        "lower_bound_count": lower_bound,
        "constructed_count": len(constructed_polynomials),
        "exact_count": exact_count,
        "matches_exact": len(constructed_polynomials) == exact_count,
    }


def _valid_s_values(n: int, k: int, ell: int) -> Iterable[int]:
    for s in range(ell, n + 1, ell):
        if s > k and s - ell < k:
            yield s


def lower_bound_grid(
    q_bits_values: Iterable[float] = (64, 128, 192, 256),
    eps_bits: float = 128.0,
    n_values: Iterable[int] = (2**10, 2**12, 2**16, 2**20),
    rates: Iterable[float] = RATES,
    mode: str = "folded",
    m: int = 1,
) -> list[dict[str, Any]]:
    """Generate numerical coset-union lower-bound rows for prize-scale parameters."""
    rows: list[dict[str, Any]] = []
    log2_budget_cache = {q_bits: q_bits - eps_bits for q_bits in q_bits_values}
    for q_bits in q_bits_values:
        for n in n_values:
            if mode == "folded" and n % m != 0:
                continue
            for rho in rates:
                k_float = rho * n
                if abs(k_float - round(k_float)) > 1e-9:
                    continue
                k = int(round(k_float))
                if not (0 < k < n):
                    continue
                threshold = threshold_params(
                    q_bits=q_bits,
                    n=n,
                    rho=rho,
                    m=m,
                    eps_bits=eps_bits,
                    mode=mode,
                )
                for ell in divisors(n):
                    for s in _valid_s_values(n, k, ell):
                        log2_list = _log2_coset_union_lower_bound(n, k, s, ell)
                        if log2_list == -math.inf:
                            continue
                        radius = 1.0 - (s / n)
                        log2_budget = log2_budget_cache[q_bits]
                        delta_entropy = float(threshold["delta_entropy"])
                        rows.append(
                            {
                                "q_bits": q_bits,
                                "n": n,
                                "k": k,
                                "rho": rho,
                                "mode": mode,
                                "m": m,
                                "ell": ell,
                                "s": s,
                                "coset_count": n // ell,
                                "cosets_chosen": s // ell,
                                "radius": radius,
                                "log2_list_lower_bound": log2_list,
                                "log2_budget": log2_budget,
                                "log2_margin": log2_list - log2_budget,
                                "beats_budget": log2_list > log2_budget,
                                "delta_entropy": delta_entropy,
                                "delta_volume_grid": float(threshold["delta_volume_grid"]),
                                "capacity_gap": delta_entropy - radius,
                                "below_capacity": radius <= delta_entropy,
                            }
                        )
    return sorted(
        rows,
        key=lambda row: (
            not bool(row["beats_budget"]),
            not bool(row["below_capacity"]),
            -float(row["radius"]),
            -float(row["log2_margin"]),
            int(row["n"]),
            float(row["rho"]),
            int(row["ell"]),
        ),
    )


LOWER_BOUND_GRID_FIELDS = [
    "q_bits",
    "n",
    "k",
    "rho",
    "mode",
    "m",
    "ell",
    "s",
    "coset_count",
    "cosets_chosen",
    "radius",
    "log2_list_lower_bound",
    "log2_budget",
    "log2_margin",
    "beats_budget",
    "delta_entropy",
    "delta_volume_grid",
    "capacity_gap",
    "below_capacity",
]


def write_lower_bound_grid_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write lower-bound grid rows to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=LOWER_BOUND_GRID_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_lower_bound_grid_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write compact Markdown summary for the lower-bound grid."""
    path.parent.mkdir(parents=True, exist_ok=True)
    beats = [row for row in rows if row["beats_budget"]]
    below = [row for row in rows if row["below_capacity"]]
    both = [row for row in rows if row["beats_budget"] and row["below_capacity"]]
    lines = [
        "# Coset-Union Lower-Bound Grid Summary",
        "",
        "Numerical scan of explicit smooth-domain coset-union lower bounds.",
        "",
        f"- Rows: `{len(rows)}`",
        f"- Rows beating `eps*q`: `{len(beats)}`",
        f"- Rows below entropy capacity candidate: `{len(below)}`",
        f"- Rows doing both: `{len(both)}`",
        "",
        "| q bits | rows | beats budget | best log2 margin | best log2 list | best radius |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for q_bits in sorted({float(row["q_bits"]) for row in rows}):
        q_rows = [row for row in rows if float(row["q_bits"]) == q_bits]
        best = max(q_rows, key=lambda row: float(row["log2_margin"]))
        lines.append(
            "| {display_q_bits:g} | {row_count} | {beats} | {log2_margin:.3f} | {log2_list_lower_bound:.3f} | {radius:.6f} |".format(
                display_q_bits=q_bits,
                row_count=len(q_rows),
                beats=sum(1 for row in q_rows if row["beats_budget"]),
                **best,
            )
        )
    lines.extend(
        [
            "",
            "| beats | below cap | radius | log2 margin | log2 list | log2 budget | q bits | n | rho | k | ell | s | delta entropy |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in rows[:40]:
        lines.append(
            "| {beats_budget} | {below_capacity} | {radius:.6f} | {log2_margin:.3f} | {log2_list_lower_bound:.3f} | {log2_budget:.3f} | {q_bits:g} | {n} | {rho:.5f} | {k} | {ell} | {s} | {delta_entropy:.6f} |".format(
                **row
            )
        )
    path.write_text("\n".join(lines) + "\n")


def _parse_number_list(text: str, cast: type = int) -> list[Any]:
    values = [part for chunk in text.split(",") for part in chunk.split()]
    return [cast(value) for value in values if value]


def _parse_rate(value: str) -> float:
    if "/" in value:
        numerator, denominator = value.split("/", 1)
        return float(numerator) / float(denominator)
    return float(value)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--q-bits", default="64,128,192,256")
    parser.add_argument("--eps-bits", type=float, default=128.0)
    parser.add_argument("--n-values", default=f"{2**10},{2**12},{2**16},{2**20}")
    parser.add_argument("--rates", default="1/2,1/4,1/8,1/16")
    parser.add_argument("--mode", default="folded", choices=("folded", "interleaved"))
    parser.add_argument("--m", type=int, default=1)
    args = parser.parse_args()

    rows = lower_bound_grid(
        q_bits_values=_parse_number_list(args.q_bits, float),
        eps_bits=args.eps_bits,
        n_values=_parse_number_list(args.n_values, int),
        rates=[_parse_rate(value) for value in _parse_number_list(args.rates, str)],
        mode=args.mode,
        m=args.m,
    )
    write_lower_bound_grid_csv(args.csv, rows)
    write_lower_bound_grid_summary(args.summary, rows)


if __name__ == "__main__":
    main()
