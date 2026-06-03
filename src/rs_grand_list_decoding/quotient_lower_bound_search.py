"""Quotient-optimized coset-union lower-bound search."""

from __future__ import annotations

import argparse
import csv
import math
import random
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable

from rs_grand_list_decoding.finite_field import factor
from rs_grand_list_decoding.lower_bound_constructions import (
    coset_union_lower_bound_count,
)
from rs_grand_list_decoding.rs_capacity_threshold import ln_comb, threshold_params


def rate_to_fraction(rho: float | str) -> tuple[int, int]:
    """Parse a rate into a reduced positive fraction."""
    if isinstance(rho, str):
        frac = Fraction(rho.strip())
    else:
        frac = Fraction(rho).limit_denominator(1_000_000)
    if frac <= 0 or frac >= 1:
        raise ValueError("rho must be in (0, 1)")
    return frac.numerator, frac.denominator


def minimal_realizing_ell(M: int, rho_num: int, rho_den: int) -> int:
    """Return the smallest ell such that rho * ell * M is integral."""
    if M <= 0:
        raise ValueError("M must be positive")
    if rho_num <= 0 or rho_den <= 0:
        raise ValueError("rate numerator and denominator must be positive")
    return rho_den // math.gcd(rho_den, M)


def quotient_candidate(M: int, rho_num: int, rho_den: int) -> dict[str, Any] | None:
    """Return quotient-level construction parameters, or None when rho*M is integral."""
    if M <= 0:
        raise ValueError("M must be positive")
    numerator = rho_num * M
    if numerator % rho_den == 0:
        return None
    r = numerator // rho_den + 1
    alpha = r / M
    return {
        "M": M,
        "rho_num": rho_num,
        "rho_den": rho_den,
        "rho": rho_num / rho_den,
        "r": r,
        "alpha": alpha,
        "radius": 1.0 - alpha,
        "log2_list_lower_bound": ln_comb(M, r) / math.log(2.0),
    }


def _max_prime_factor(n: int) -> int:
    if n <= 1:
        return 1
    return max(factor(n))


def is_B_smooth(n: int, B: int) -> bool:
    """Return whether every prime factor of n is at most B."""
    if n <= 0:
        raise ValueError("n must be positive")
    if B < 2:
        return n == 1
    return _max_prime_factor(n) <= B


def _threshold_mode(mode: str) -> str:
    if mode == "scalar":
        return "folded"
    if mode in {"folded", "interleaved"}:
        return mode
    raise ValueError("mode must be 'scalar', 'folded', or 'interleaved'")


def generate_quotient_candidates(
    q_bits_values: Iterable[float] = (128, 192, 256),
    eps_bits: float = 128.0,
    rates: Iterable[str | float] = ("1/2", "1/4", "1/8", "1/16"),
    M_max: int = 5000,
    smoothness_bounds: Iterable[int] = (16, 32, 64, 128),
    n_scale_multipliers: Iterable[int] = (1, 2, 4, 8, 16, 32, 64),
    mode: str = "scalar",
) -> list[dict[str, Any]]:
    """Generate smooth quotient-level coset-union lower-bound candidates."""
    smoothness_bounds = tuple(sorted(set(int(bound) for bound in smoothness_bounds)))
    if not smoothness_bounds:
        raise ValueError("at least one smoothness bound is required")
    max_smoothness = max(smoothness_bounds)
    multipliers = tuple(sorted(set(int(multiplier) for multiplier in n_scale_multipliers)))
    if any(multiplier <= 0 for multiplier in multipliers):
        raise ValueError("n scale multipliers must be positive")

    parsed_rates = [rate_to_fraction(rho) for rho in rates]
    threshold_mode = _threshold_mode(mode)
    threshold_cache: dict[tuple[float, int, float], dict[str, float | int | str]] = {}
    rows: list[dict[str, Any]] = []

    for rho_num, rho_den in parsed_rates:
        rho = rho_num / rho_den
        for M in range(1, M_max + 1):
            candidate = quotient_candidate(M, rho_num, rho_den)
            if candidate is None:
                continue
            M_smooth_B = _max_prime_factor(M)
            if M_smooth_B > max_smoothness:
                continue
            ell0 = minimal_realizing_ell(M, rho_num, rho_den)
            r = int(candidate["r"])
            radius = float(candidate["radius"])
            log2_list = float(candidate["log2_list_lower_bound"])
            for multiplier in multipliers:
                ell = ell0 * multiplier
                n = ell * M
                n_smooth_B = _max_prime_factor(n)
                if n_smooth_B > max_smoothness:
                    continue
                k = (rho_num * n) // rho_den
                if rho_num * n % rho_den != 0:
                    raise AssertionError("minimal ell did not realize integral k")
                s = r * ell
                if coset_union_lower_bound_count(n, k, s, ell) == 0:
                    raise AssertionError("quotient candidate failed theorem conditions")
                for q_bits in q_bits_values:
                    q_bits = float(q_bits)
                    threshold_key = (q_bits, n, rho)
                    if threshold_key not in threshold_cache:
                        threshold_cache[threshold_key] = threshold_params(
                            q_bits=q_bits,
                            n=n,
                            rho=rho,
                            m=1,
                            eps_bits=eps_bits,
                            mode=threshold_mode,
                        )
                    threshold = threshold_cache[threshold_key]
                    delta_entropy = float(threshold["delta_entropy"])
                    log2_budget = q_bits - eps_bits
                    requested_smooth_B = next(
                        bound for bound in smoothness_bounds if n_smooth_B <= bound
                    )
                    rows.append(
                        {
                            "q_bits": q_bits,
                            "eps_bits": eps_bits,
                            "rho": rho,
                            "rho_num": rho_num,
                            "rho_den": rho_den,
                            "M": M,
                            "r": r,
                            "alpha": float(candidate["alpha"]),
                            "radius": radius,
                            "ell0": ell0,
                            "ell_multiplier": multiplier,
                            "ell": ell,
                            "n": n,
                            "k": k,
                            "s": s,
                            "mode": mode,
                            "M_smooth_B": M_smooth_B,
                            "n_smooth_B": n_smooth_B,
                            "requested_smooth_B": requested_smooth_B,
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
            -float(row["log2_margin"]),
            -float(row["radius"]),
            float(row["q_bits"]),
            float(row["rho"]),
            int(row["M"]),
            int(row["ell"]),
        ),
    )


def _decompose_odd_part(n: int) -> tuple[int, int]:
    d = n - 1
    s = 0
    while d % 2 == 0:
        s += 1
        d //= 2
    return s, d


def is_probable_prime(n: int, rounds: int = 32) -> bool:
    """Miller-Rabin primality test, deterministic below 2^64."""
    if n < 2:
        return False
    small_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    if n in small_primes:
        return True
    if any(n % prime == 0 for prime in small_primes):
        return False

    if n < 2**64:
        bases = (2, 325, 9375, 28178, 450775, 9_780_504, 1_795_265_022)
    else:
        rng = random.Random(n)
        bases = small_primes + tuple(rng.randrange(2, n - 2) for _ in range(rounds))

    s, d = _decompose_odd_part(n)
    for a in bases:
        a %= n
        if a in (0, 1):
            continue
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True


def find_prime_for_domain(n: int, q_bits: int, max_trials: int = 100_000) -> int | None:
    """Search for a prime p = c*n + 1 with bit length about q_bits."""
    if n <= 0:
        raise ValueError("n must be positive")
    if q_bits <= 1:
        raise ValueError("q_bits must exceed 1")
    lower = 1 << (q_bits - 1)
    c = max(1, (lower - 1 + n - 1) // n)
    for _ in range(max_trials):
        p = c * n + 1
        if p.bit_length() > q_bits:
            return None
        if p.bit_length() == q_bits and is_probable_prime(p):
            return p
        c += 1
    return None


def materialize_candidate(row: dict[str, Any], max_trials: int = 100_000) -> dict[str, Any]:
    """Search for a concrete prime field realizing a quotient candidate row."""
    q_bits = int(float(row["q_bits"]))
    n = int(row["n"])
    p = find_prime_for_domain(n, q_bits=q_bits, max_trials=max_trials)
    out = {
        "q_bits": q_bits,
        "rho": row["rho"],
        "M": row["M"],
        "r": row["r"],
        "ell": row["ell"],
        "n": n,
        "k": row["k"],
        "s": row["s"],
        "center": f"x^{row['s']}",
        "lower_bound": f"binom({row['M']},{row['r']})",
        "radius": row["radius"],
        "delta_entropy": row["delta_entropy"],
        "log2_list_lower_bound": row["log2_list_lower_bound"],
        "log2_budget": row["log2_budget"],
        "log2_margin": row["log2_margin"],
        "M_smooth_B": row["M_smooth_B"],
        "n_smooth_B": row["n_smooth_B"],
        "below_capacity": row["below_capacity"],
        "beats_budget": row["beats_budget"],
        "prime_found": p is not None,
        "p": p or "",
    }
    return out


def select_materialization_rows(rows: list[dict[str, Any]], limit: int) -> list[dict[str, Any]]:
    """Select top rows plus canonical Stage 7 examples for prime materialization."""
    selected: list[dict[str, Any]] = []
    seen: set[tuple[Any, ...]] = set()

    def add(row: dict[str, Any]) -> None:
        key = (row["q_bits"], row["rho_num"], row["rho_den"], row["M"], row["ell"])
        if key in seen:
            return
        seen.add(key)
        selected.append(row)

    for row in rows[:limit]:
        add(row)
    for row in rows:
        if (
            int(float(row["q_bits"])) == 256
            and int(row["rho_num"]) == 1
            and int(row["rho_den"]) == 4
            and int(row["M"]) == 225
            and int(row["ell_multiplier"]) == 1
        ):
            add(row)
            break
    return selected


QUOTIENT_FIELDS = [
    "q_bits",
    "eps_bits",
    "rho",
    "rho_num",
    "rho_den",
    "M",
    "r",
    "alpha",
    "radius",
    "ell0",
    "ell_multiplier",
    "ell",
    "n",
    "k",
    "s",
    "mode",
    "M_smooth_B",
    "n_smooth_B",
    "requested_smooth_B",
    "log2_list_lower_bound",
    "log2_budget",
    "log2_margin",
    "beats_budget",
    "delta_entropy",
    "delta_volume_grid",
    "capacity_gap",
    "below_capacity",
]


MATERIALIZED_FIELDS = [
    "q_bits",
    "rho",
    "M",
    "r",
    "ell",
    "n",
    "k",
    "s",
    "center",
    "lower_bound",
    "radius",
    "delta_entropy",
    "log2_list_lower_bound",
    "log2_budget",
    "log2_margin",
    "M_smooth_B",
    "n_smooth_B",
    "below_capacity",
    "beats_budget",
    "prime_found",
    "p",
]


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    """Write rows to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, rows: list[dict[str, Any]], materialized: list[dict[str, Any]]) -> None:
    """Write a compact Markdown summary of quotient candidates."""
    path.parent.mkdir(parents=True, exist_ok=True)
    beats = [row for row in rows if row["beats_budget"]]
    below = [row for row in rows if row["below_capacity"]]
    both = [row for row in rows if row["beats_budget"] and row["below_capacity"]]
    q256_both = [row for row in both if int(float(row["q_bits"])) == 256]
    lines = [
        "# Quotient Lower-Bound Search Summary",
        "",
        "Quotient-optimized scan of the coset-union lower-bound construction over mixed-smooth domains.",
        "",
        f"- Rows: `{len(rows)}`",
        f"- Rows beating `eps*q`: `{len(beats)}`",
        f"- Rows below entropy capacity candidate: `{len(below)}`",
        f"- Rows doing both: `{len(both)}`",
        f"- `q_bits=256` rows doing both: `{len(q256_both)}`",
        "",
        "| margin | radius | delta entropy | q bits | rho | M | r | ell | n | k | s | M smooth | n smooth |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows[:40]:
        lines.append(
            "| {log2_margin:.3f} | {radius:.6f} | {delta_entropy:.6f} | {q_bits:g} | {rho:.5f} | {M} | {r} | {ell} | {n} | {k} | {s} | {M_smooth_B} | {n_smooth_B} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Materialized Primes",
            "",
            "| found | q bits | rho | M | ell | n | log2 margin | p |",
            "|---|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in materialized:
        p_text = str(row["p"])
        if len(p_text) > 36:
            p_text = p_text[:18] + "..." + p_text[-12:]
        lines.append(
            "| {prime_found} | {q_bits} | {rho:.5f} | {M} | {ell} | {n} | {log2_margin:.3f} | `{p}` |".format(
                **{**row, "p": p_text}
            )
        )
    path.write_text("\n".join(lines) + "\n")


def _parse_number_list(text: str, cast: type = int) -> list[Any]:
    values = [part for chunk in text.split(",") for part in chunk.split()]
    return [cast(value) for value in values if value]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--materialized-csv", type=Path, required=True)
    parser.add_argument("--q-bits", default="128,192,256")
    parser.add_argument("--eps-bits", type=float, default=128.0)
    parser.add_argument("--rates", default="1/2,1/4,1/8,1/16")
    parser.add_argument("--M-max", type=int, default=5000)
    parser.add_argument("--smoothness-bounds", default="16,32,64,128")
    parser.add_argument("--n-scale-multipliers", default="1,2,4,8,16,32,64")
    parser.add_argument("--mode", default="scalar")
    parser.add_argument("--materialize-limit", type=int, default=20)
    parser.add_argument("--prime-max-trials", type=int, default=100_000)
    args = parser.parse_args()

    rows = generate_quotient_candidates(
        q_bits_values=_parse_number_list(args.q_bits, float),
        eps_bits=args.eps_bits,
        rates=_parse_number_list(args.rates, str),
        M_max=args.M_max,
        smoothness_bounds=_parse_number_list(args.smoothness_bounds, int),
        n_scale_multipliers=_parse_number_list(args.n_scale_multipliers, int),
        mode=args.mode,
    )
    materialization_rows = select_materialization_rows(rows, args.materialize_limit)
    materialized = [
        materialize_candidate(row, max_trials=args.prime_max_trials)
        for row in materialization_rows
    ]
    write_csv(args.csv, rows, QUOTIENT_FIELDS)
    write_csv(args.materialized_csv, materialized, MATERIALIZED_FIELDS)
    write_summary(args.summary, rows, materialized)


if __name__ == "__main__":
    main()
