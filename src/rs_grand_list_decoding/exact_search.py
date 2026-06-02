"""Exact tiny-parameter max-list search for Reed-Solomon experiments."""

from __future__ import annotations

import argparse
import csv
import itertools
import math
from pathlib import Path
from typing import Any, Iterable

from rs_grand_list_decoding.finite_field import factor, smooth_domain
from rs_grand_list_decoding.rs_capacity_threshold import threshold_params
from rs_grand_list_decoding.rs_code import (
    enumerate_rs_codewords,
    fold_codeword,
    hamming_distance,
    interleave_codewords,
)


DEFAULT_MAX_WORDS = 2_000_000
DEFAULT_MAX_CODEWORDS = 200_000
DEFAULT_MAX_OPERATIONS = 5_000_000


def alphabet_from_codewords(codewords: list[tuple]) -> list[Any]:
    """Return sorted alphabet symbols appearing in codewords."""
    symbols = {symbol for word in codewords for symbol in word}
    return sorted(symbols, key=repr)


def received_space_size(alphabet: list[Any], N: int) -> int:
    """Return |alphabet|^N."""
    if N < 0:
        raise ValueError("N must be nonnegative")
    return len(alphabet) ** N


def full_tuple_space(
    alphabet: list[Any],
    N: int,
    max_words: int = DEFAULT_MAX_WORDS,
    override: bool = False,
) -> Iterable[tuple]:
    """Enumerate all received words in alphabet^N, refusing infeasible cases by default."""
    total = received_space_size(alphabet, N)
    if total > max_words and not override:
        raise ValueError(
            f"received-word space has {total} words, exceeding max_words={max_words}"
        )
    return itertools.product(alphabet, repeat=N)


def _check_search_budget(
    alphabet: list[Any],
    N: int,
    codeword_count: int,
    max_words: int,
    max_operations: int,
    override: bool,
) -> None:
    total_words = received_space_size(alphabet, N)
    if total_words > max_words and not override:
        raise ValueError(
            f"received-word space has {total_words} words, exceeding max_words={max_words}"
        )
    operations = total_words * codeword_count
    if operations > max_operations and not override:
        raise ValueError(
            f"search would require {operations} center/codeword distance checks, "
            f"exceeding max_operations={max_operations}"
        )


def _validate_codewords(codewords: list[tuple]) -> int:
    if not codewords:
        raise ValueError("codewords must be nonempty")
    N = len(codewords[0])
    if any(len(word) != N for word in codewords):
        raise ValueError("all codewords must have the same length")
    return N


def exact_max_list(
    codewords: list[tuple],
    alphabet: list[Any],
    radius: int,
    max_words: int = DEFAULT_MAX_WORDS,
    max_operations: int = DEFAULT_MAX_OPERATIONS,
    override: bool = False,
) -> dict[str, Any]:
    """Compute the exact maximum list size at one integer radius."""
    N = _validate_codewords(codewords)
    if radius < 0 or radius > N:
        raise ValueError("radius must be in [0, N]")
    _check_search_budget(
        alphabet,
        N,
        len(codewords),
        max_words=max_words,
        max_operations=max_operations,
        override=override,
    )

    max_count = -1
    best_center: tuple | None = None
    num_centers = 0
    for center in full_tuple_space(alphabet, N, max_words=max_words, override=override):
        count = sum(hamming_distance(center, codeword) <= radius for codeword in codewords)
        if count > max_count:
            max_count = count
            best_center = center
            num_centers = 1
        elif count == max_count:
            num_centers += 1

    return {
        "radius": radius,
        "max_list": max_count,
        "center": best_center,
        "num_centers": num_centers,
    }


def exact_profile(
    codewords: list[tuple],
    alphabet: list[Any],
    max_words: int = DEFAULT_MAX_WORDS,
    max_operations: int = DEFAULT_MAX_OPERATIONS,
    override: bool = False,
) -> list[dict[str, Any]]:
    """Compute exact maximum list size for every radius 0..N."""
    N = _validate_codewords(codewords)
    _check_search_budget(
        alphabet,
        N,
        len(codewords),
        max_words=max_words,
        max_operations=max_operations,
        override=override,
    )
    profile = [
        {"radius": radius, "max_list": -1, "center": None, "num_centers": 0}
        for radius in range(N + 1)
    ]

    for center in full_tuple_space(alphabet, N, max_words=max_words, override=override):
        distances = sorted(hamming_distance(center, codeword) for codeword in codewords)
        count = 0
        index = 0
        for radius in range(N + 1):
            while index < len(distances) and distances[index] <= radius:
                count += 1
                index += 1
            row = profile[radius]
            if count > row["max_list"]:
                row["max_list"] = count
                row["center"] = center
                row["num_centers"] = 1
            elif count == row["max_list"]:
                row["num_centers"] += 1
    return profile


def build_scalar_rs_instance(p: int, n: int, k: int) -> dict[str, Any]:
    """Build a tiny scalar RS instance over a smooth subgroup domain."""
    domain = smooth_domain(p, n)
    codewords = enumerate_rs_codewords(p, domain, k)
    return {
        "mode": "scalar",
        "p": p,
        "n": n,
        "k": k,
        "m": 1,
        "domain": domain,
        "codewords": codewords,
        "alphabet": list(range(p)),
        "ambient_alphabet_size": p,
    }


def build_folded_rs_instance(
    p: int,
    n: int,
    k: int,
    m: int,
    full_alphabet: bool = False,
) -> dict[str, Any]:
    """Build a folded scalar RS instance."""
    scalar = build_scalar_rs_instance(p, n, k)
    codewords = [fold_codeword(word, m) for word in scalar["codewords"]]
    alphabet = (
        list(itertools.product(range(p), repeat=m))
        if full_alphabet
        else alphabet_from_codewords(codewords)
    )
    return {
        "mode": "folded",
        "p": p,
        "n": n,
        "k": k,
        "m": m,
        "domain": scalar["domain"],
        "codewords": codewords,
        "alphabet": alphabet,
        "ambient_alphabet_size": p**m,
    }


def build_interleaved_rs_instance(
    p: int,
    n: int,
    k: int,
    m: int,
    full_alphabet: bool = False,
    max_codewords: int = DEFAULT_MAX_CODEWORDS,
    override: bool = False,
) -> dict[str, Any]:
    """Build a direct interleaved RS instance by enumerating m-tuples of scalar codewords."""
    scalar = build_scalar_rs_instance(p, n, k)
    scalar_codewords = scalar["codewords"]
    total_codewords = len(scalar_codewords) ** m
    if total_codewords > max_codewords and not override:
        raise ValueError(
            f"interleaved code has {total_codewords} codewords, "
            f"exceeding max_codewords={max_codewords}"
        )
    codewords = [
        interleave_codewords(list(words))
        for words in itertools.product(scalar_codewords, repeat=m)
    ]
    alphabet = (
        list(itertools.product(range(p), repeat=m))
        if full_alphabet
        else alphabet_from_codewords(codewords)
    )
    return {
        "mode": "interleaved",
        "p": p,
        "n": n,
        "k": k,
        "m": m,
        "domain": scalar["domain"],
        "codewords": codewords,
        "alphabet": alphabet,
        "ambient_alphabet_size": p**m,
    }


def divisors(n: int) -> list[int]:
    """Return positive divisors of n."""
    factors = factor(n)
    values = [1]
    for prime, exponent in factors.items():
        values = [
            value * prime_power
            for value in values
            for prime_power in (prime**power for power in range(exponent + 1))
        ]
    return sorted(values)


def symbol_min_distance(mode: str, n: int, k: int, m: int) -> int:
    """Return the symbol-level minimum distance for scalar/interleaved/folded modes."""
    if mode == "folded":
        if n % m != 0:
            raise ValueError("folded mode requires m | n")
        return (n // m) - ((k - 1) // m)
    return n - k + 1


def _center_to_string(center: tuple | None, max_len: int = 120) -> str:
    if center is None:
        return ""
    text = repr(center)
    if len(text) <= max_len:
        return text
    return text[: max_len - 3] + "..."


def _thresholds_for_instance(instance: dict[str, Any], eps_bits: float) -> dict[str, Any]:
    mode = instance["mode"]
    if mode == "scalar":
        threshold_mode = "folded"
        m = 1
    else:
        threshold_mode = mode
        m = instance["m"]
    out = threshold_params(
        q_bits=math.log2(instance["p"]),
        n=instance["n"],
        rho=instance["k"] / instance["n"],
        m=m,
        eps_bits=eps_bits,
        mode=threshold_mode,
    )
    return {
        "volume_grid": out["delta_volume_grid"],
        "capacity_entropy": out["delta_entropy"],
    }


def _instance_rows(
    instance: dict[str, Any],
    eps_bits: float,
    max_words: int,
    max_operations: int,
) -> list[dict[str, Any]]:
    N = len(instance["codewords"][0])
    alphabet = instance["alphabet"]
    rho = instance["k"] / instance["n"]
    base = {
        "p": instance["p"],
        "n": instance["n"],
        "k": instance["k"],
        "rho": rho,
        "m": instance["m"],
        "mode": instance["mode"],
        "N": N,
        "search_alphabet_size": len(alphabet),
        "ambient_alphabet_size": instance["ambient_alphabet_size"],
        "code_size": len(instance["codewords"]),
        "johnson_radius": 1.0 - math.sqrt(rho),
        **_thresholds_for_instance(instance, eps_bits),
    }

    try:
        profile = exact_profile(
            instance["codewords"],
            alphabet,
            max_words=max_words,
            max_operations=max_operations,
        )
    except ValueError as exc:
        return [
            {
                **base,
                "radius": "",
                "max_list": "",
                "num_centers": "",
                "best_center": "",
                "status": f"refused: {exc}",
            }
        ]

    return [
        {
            **base,
            "radius": row["radius"],
            "max_list": row["max_list"],
            "num_centers": row["num_centers"],
            "best_center": _center_to_string(row["center"]),
            "status": "ok",
        }
        for row in profile
    ]


def generate_exact_sweep_rows(
    primes: Iterable[int] = (5, 7, 13, 17),
    max_n: int = 8,
    m_values: Iterable[int] = (1, 2),
    eps_bits: float = 0.0,
    max_words: int = DEFAULT_MAX_WORDS,
    max_operations: int = DEFAULT_MAX_OPERATIONS,
) -> list[dict[str, Any]]:
    """Generate the first exact tiny-parameter sweep requested in the spec."""
    rows: list[dict[str, Any]] = []
    for p in primes:
        for n in divisors(p - 1):
            if n <= 1 or n > max_n:
                continue
            for k in range(1, n):
                if p**k > DEFAULT_MAX_CODEWORDS:
                    rho = k / n
                    rows.append(
                        {
                            "p": p,
                            "n": n,
                            "k": k,
                            "rho": rho,
                            "m": 1,
                            "mode": "scalar",
                            "N": n,
                            "search_alphabet_size": p,
                            "ambient_alphabet_size": p,
                            "code_size": "",
                            "radius": "",
                            "max_list": "",
                            "num_centers": "",
                            "johnson_radius": 1.0 - math.sqrt(rho),
                            "volume_grid": "",
                            "capacity_entropy": "",
                            "best_center": "",
                            "status": (
                                f"refused: scalar code has {p**k} codewords, "
                                f"exceeding max_codewords={DEFAULT_MAX_CODEWORDS}"
                            ),
                        }
                    )
                    continue
                scalar = build_scalar_rs_instance(p, n, k)
                rows.extend(
                    _instance_rows(
                        scalar,
                        eps_bits=eps_bits,
                        max_words=max_words,
                        max_operations=max_operations,
                    )
                )

                for m in m_values:
                    if m <= 1:
                        continue
                    if n % m == 0:
                        folded = build_folded_rs_instance(p, n, k, m)
                        rows.extend(
                            _instance_rows(
                                folded,
                                eps_bits=eps_bits,
                                max_words=max_words,
                                max_operations=max_operations,
                            )
                        )
                    try:
                        interleaved = build_interleaved_rs_instance(p, n, k, m)
                    except ValueError as exc:
                        N = n
                        rho = k / n
                        rows.append(
                            {
                                "p": p,
                                "n": n,
                                "k": k,
                                "rho": rho,
                                "m": m,
                                "mode": "interleaved",
                                "N": N,
                                "search_alphabet_size": "",
                                "ambient_alphabet_size": p**m,
                                "code_size": "",
                                "radius": "",
                                "max_list": "",
                                "num_centers": "",
                                "johnson_radius": 1.0 - math.sqrt(rho),
                                "volume_grid": "",
                                "capacity_entropy": "",
                                "best_center": "",
                                "status": f"refused: {exc}",
                            }
                        )
                    else:
                        rows.extend(
                            _instance_rows(
                                interleaved,
                                eps_bits=eps_bits,
                                max_words=max_words,
                                max_operations=max_operations,
                            )
                        )
    return rows


def write_exact_sweep_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write exact-sweep rows to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "p",
        "n",
        "k",
        "rho",
        "m",
        "mode",
        "N",
        "search_alphabet_size",
        "ambient_alphabet_size",
        "code_size",
        "radius",
        "max_list",
        "num_centers",
        "johnson_radius",
        "volume_grid",
        "capacity_entropy",
        "best_center",
        "status",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_exact_sweep_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write a compact Markdown summary of the exact sweep."""
    path.parent.mkdir(parents=True, exist_ok=True)
    ok_rows = [row for row in rows if row["status"] == "ok"]
    refused_rows = [row for row in rows if str(row["status"]).startswith("refused")]

    grouped: dict[tuple, list[dict[str, Any]]] = {}
    for row in ok_rows:
        key = (row["p"], row["n"], row["k"], row["mode"], row["m"])
        grouped.setdefault(key, []).append(row)

    interesting_rows: list[tuple[str, dict[str, Any]]] = []
    for group_rows in grouped.values():
        first = group_rows[0]
        N = int(first["N"])
        n = int(first["n"])
        k = int(first["k"])
        m = int(first["m"])
        mode = str(first["mode"])
        rho = float(first["rho"])
        d = symbol_min_distance(mode, n, k, m)
        candidate_radii = {
            "johnson": math.floor((1.0 - math.sqrt(rho)) * N),
            "capacity": math.floor(float(first["capacity_entropy"]) * N),
            "d-1": max(0, min(N, d - 1)),
            "d": max(0, min(N, d)),
        }
        by_radius = {row["radius"]: row for row in group_rows}
        for label, radius in candidate_radii.items():
            if radius in by_radius:
                interesting_rows.append((label, by_radius[radius]))

    lines = [
        "# Exact Sweep Summary",
        "",
        "Toy exact enumeration over smooth subgroup domains.",
        "",
        f"- Successful profile rows: `{len(ok_rows)}`",
        f"- Refused rows: `{len(refused_rows)}`",
        "- Default sweep includes `p in {5, 7, 13, 17}` and smooth divisors `n <= 8`.",
        "- Toy capacity columns use `eps_bits = 0`, not the prize value `2^-128`.",
        "",
        "| label | p | n | k | mode | m | radius | max list | search Q | ambient Q | code size |",
        "|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|",
    ]
    for label, row in interesting_rows[:60]:
        lines.append(
            "| {label} | {p} | {n} | {k} | {mode} | {m} | {radius} | {max_list} | {search_alphabet_size} | {ambient_alphabet_size} | {code_size} |".format(
                label=label, **row
            )
        )
    if len(interesting_rows) > 60:
        lines.append(
            "| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | "
            f"{len(interesting_rows) - 60} more rows |"
        )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--max-words", type=int, default=DEFAULT_MAX_WORDS)
    parser.add_argument("--max-operations", type=int, default=DEFAULT_MAX_OPERATIONS)
    parser.add_argument("--eps-bits", type=float, default=0.0)
    args = parser.parse_args()

    rows = generate_exact_sweep_rows(
        max_words=args.max_words,
        max_operations=args.max_operations,
        eps_bits=args.eps_bits,
    )
    write_exact_sweep_csv(args.csv, rows)
    write_exact_sweep_summary(args.summary, rows)


if __name__ == "__main__":
    main()
