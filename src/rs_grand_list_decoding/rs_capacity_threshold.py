"""Capacity and volume thresholds for the RS grand list-decoding challenge."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


RATES = (1 / 2, 1 / 4, 1 / 8, 1 / 16)


def _ln_q_minus_1(ln_q: float) -> float:
    """Compute ln(Q - 1) from ln(Q) without overflowing for huge Q."""
    if ln_q > 50:
        return ln_q
    return math.log(math.expm1(ln_q))


def h_q(delta: float, ln_q: float) -> float:
    """Return the Q-ary entropy H_Q(delta), where ln_q is ln(Q)."""
    if ln_q <= 0:
        raise ValueError("ln_q must be positive")
    if delta <= 0:
        return 0.0
    if delta >= 1:
        return _ln_q_minus_1(ln_q) / ln_q

    binary_entropy_nats = -(
        delta * math.log(delta) + (1.0 - delta) * math.log(1.0 - delta)
    )
    return delta * (_ln_q_minus_1(ln_q) / ln_q) + binary_entropy_nats / ln_q


def inverse_h_q(target: float, ln_q: float, iters: int = 100) -> float:
    """Invert H_Q on its increasing branch by bisection."""
    if ln_q <= 0:
        raise ValueError("ln_q must be positive")
    if target <= 0:
        return 0.0

    hi = 1.0 if ln_q > 50 else 1.0 - math.exp(-ln_q)
    max_entropy = h_q(hi, ln_q)
    if target >= max_entropy:
        return hi

    lo = 0.0
    for _ in range(iters):
        mid = (lo + hi) / 2.0
        if h_q(mid, ln_q) <= target:
            lo = mid
        else:
            hi = mid
    return lo


def log_add_exp(a: float, b: float) -> float:
    """Stable log(exp(a) + exp(b))."""
    if a == -math.inf:
        return b
    if b == -math.inf:
        return a
    if a < b:
        a, b = b, a
    return a + math.log1p(math.exp(b - a))


def ln_comb(n: int, k: int) -> float:
    """Return ln(binomial(n, k))."""
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def log_q_ball(N: int, t: int, ln_q: float) -> float:
    """Return log_Q V_Q(N, t), computed in natural-log space."""
    if N < 0:
        raise ValueError("N must be nonnegative")
    if t < 0:
        return -math.inf
    if t > N:
        t = N
    if ln_q <= 0:
        raise ValueError("ln_q must be positive")
    if t == N:
        return float(N)

    ln_q_minus_1 = _ln_q_minus_1(ln_q)

    # For large alphabets, the outer shell dominates the Hamming ball. Summing
    # backwards from the outer shell avoids O(N) work for SNARK-sized parameters.
    if ln_q >= 20:
        ln_top = ln_comb(N, t) + t * ln_q_minus_1
        log_relative_total = 0.0
        log_relative_term = 0.0
        for i in range(t, 0, -1):
            log_relative_term += (
                math.log(i) - math.log(N - i + 1) - ln_q_minus_1
            )
            log_relative_total = log_add_exp(log_relative_total, log_relative_term)
            if log_relative_term < -60:
                break
        return (ln_top + log_relative_total) / ln_q

    total = -math.inf
    for i in range(t + 1):
        term = ln_comb(N, i) + i * ln_q_minus_1
        total = log_add_exp(total, term)
    return total / ln_q


@dataclass(frozen=True)
class ThresholdParams:
    """Computed threshold parameters for one configuration."""

    mode: str
    q_bits: float
    n: int
    m: int
    rho: float
    N: int
    Q_bits: float
    rate: float
    entropy_target: float
    delta_entropy: float
    volume_grid_t: int
    delta_volume_grid: float

    def as_dict(self) -> dict[str, float | int | str]:
        return {
            "mode": self.mode,
            "q_bits": self.q_bits,
            "n": self.n,
            "m": self.m,
            "rho": self.rho,
            "N": self.N,
            "Q_bits": self.Q_bits,
            "rate": self.rate,
            "entropy_target": self.entropy_target,
            "delta_entropy": self.delta_entropy,
            "volume_grid_t": self.volume_grid_t,
            "delta_volume_grid": self.delta_volume_grid,
        }


def threshold_params(
    q_bits: float,
    n: int,
    rho: float,
    m: int,
    eps_bits: float = 128.0,
    mode: str = "folded",
) -> dict[str, float | int | str]:
    """
    Compute entropy and exact-grid volume thresholds.

    `q = 2^q_bits` and `eps = 2^-eps_bits`.

    Modes:
    - folded: N = n / m, Q = q^m, R = rho
    - interleaved: N = n, Q = q^m, R = rho
    """
    if q_bits <= 0:
        raise ValueError("q_bits must be positive")
    if n <= 0:
        raise ValueError("n must be positive")
    if m <= 0:
        raise ValueError("m must be positive")
    if not (0.0 < rho < 1.0):
        raise ValueError("rho must be in (0, 1)")

    if mode == "folded":
        if n % m != 0:
            raise ValueError("folded mode requires m | n")
        N = n // m
        rate = rho
    elif mode == "interleaved":
        N = n
        rate = rho
    else:
        raise ValueError("mode must be 'folded' or 'interleaved'")

    ln_q = q_bits * math.log(2.0)
    ln_Q = m * ln_q
    Q_bits = m * q_bits

    logQ_eps_q = (q_bits - eps_bits) / Q_bits
    entropy_target = 1.0 - rate + logQ_eps_q / N
    delta_entropy = inverse_h_q(entropy_target, ln_Q)

    rhs = N * (1.0 - rate) + logQ_eps_q
    lo, hi = 0, N
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if log_q_ball(N, mid, ln_Q) <= rhs:
            lo = mid
        else:
            hi = mid - 1

    result = ThresholdParams(
        mode=mode,
        q_bits=q_bits,
        n=n,
        m=m,
        rho=rho,
        N=N,
        Q_bits=Q_bits,
        rate=rate,
        entropy_target=entropy_target,
        delta_entropy=delta_entropy,
        volume_grid_t=lo,
        delta_volume_grid=lo / N,
    )
    return result.as_dict()


def generate_capacity_rows(
    q_bits_values: Iterable[float] = (64, 128, 192, 256),
    n_values: Iterable[int] = (2**10, 2**16, 2**20, 2**24),
    m_values: Iterable[int] = (1, 2, 4, 8),
    rates: Iterable[float] = RATES,
    modes: Iterable[str] = ("folded", "interleaved"),
    eps_bits: float = 128.0,
) -> list[dict[str, float | int | str]]:
    """Generate the capacity-table experiment rows from SPEC.md."""
    rows: list[dict[str, float | int | str]] = []
    for q_bits in q_bits_values:
        for n in n_values:
            for m in m_values:
                for rho in rates:
                    for mode in modes:
                        if mode == "folded" and n % m != 0:
                            continue
                        rows.append(
                            threshold_params(
                                q_bits=q_bits,
                                n=n,
                                rho=rho,
                                m=m,
                                eps_bits=eps_bits,
                                mode=mode,
                            )
                        )
    return rows


def write_capacity_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    """Write capacity rows to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "q_bits",
        "n",
        "m",
        "rho",
        "mode",
        "N",
        "Q_bits",
        "rate",
        "entropy_target",
        "delta_entropy",
        "volume_grid_t",
        "delta_volume_grid",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    """Write a compact Markdown summary for the most cited q=2^256 cases."""
    path.parent.mkdir(parents=True, exist_ok=True)
    selected = [
        row
        for row in rows
        if row["q_bits"] == 256
        and row["n"] == 2**20
        and row["mode"] == "folded"
        and row["m"] in (1, 4)
    ]
    selected.sort(key=lambda r: (int(r["m"]), -float(r["rho"])))

    lines = [
        "# Capacity Table Summary",
        "",
        "Parameters: `q = 2^256`, `n = 2^20`, `eps = 2^-128`, folded mode.",
        "",
        "| m | rho | delta_entropy | delta_volume_grid | volume_grid_t |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in selected:
        lines.append(
            "| {m} | {rho:.5f} | {delta_entropy:.8f} | {delta_volume_grid:.8f} | {volume_grid_t} |".format(
                m=int(row["m"]),
                rho=float(row["rho"]),
                delta_entropy=float(row["delta_entropy"]),
                delta_volume_grid=float(row["delta_volume_grid"]),
                volume_grid_t=int(row["volume_grid_t"]),
            )
        )
    lines.append("")
    lines.append(
        "The exact grid value is the largest integer radius satisfying the universal volume converse."
    )
    path.write_text("\n".join(lines) + "\n")


def _print_default_table() -> None:
    print("rho,delta_entropy,delta_volume_grid")
    for rho in RATES:
        out = threshold_params(q_bits=256, n=2**20, rho=rho, m=1, mode="folded")
        print(f"{rho:.5f},{out['delta_entropy']:.8f},{out['delta_volume_grid']:.8f}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=Path, help="write the full capacity table to this CSV")
    parser.add_argument("--summary", type=Path, help="write a Markdown summary")
    args = parser.parse_args()

    if args.csv or args.summary:
        rows = generate_capacity_rows()
        if args.csv:
            write_capacity_csv(args.csv, rows)
        if args.summary:
            write_summary(args.summary, rows)
    else:
        _print_default_table()


if __name__ == "__main__":
    main()
