#!/usr/bin/env python3
"""Plot Design 2 risk and CI summaries from sim2.R CSV outputs."""

from __future__ import annotations

import csv
import os
import re
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="matplotlib-"))

import matplotlib.pyplot as plt


FILE_RE = re.compile(
    r"^(risk|CI_len|CI_coverage)_(\d+)_([0-9.]+)_([0-9.]+)_(\d+)_(\d+)\.csv$"
)
METHODS = ["MMA", "AIC", "BIC", "SAIC", "SBIC", "FULL", "LASSO"]
COVERAGE_METHODS = METHODS + ["unadj"]
COLORS = {
    "MMA": "black",
    "AIC": "tab:blue",
    "BIC": "tab:green",
    "SAIC": "tab:orange",
    "SBIC": "tab:red",
    "FULL": "brown",
    "LASSO": "tab:purple",
    "unadj": "0.45",
}
MARKERS = {
    "MMA": "s",
    "AIC": "^",
    "BIC": "x",
    "SAIC": "v",
    "SBIC": "D",
    "FULL": "o",
    "LASSO": "*",
    "unadj": "+",
}


def read_summary(path: Path) -> tuple[list[float], dict[str, list[float]]]:
    with path.open(newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        methods = header[1:]
        x_values: list[float] = []
        data = {method: [] for method in methods}
        for row in reader:
            if not row:
                continue
            x_values.append(float(row[0]))
            for method, value in zip(methods, row[1:]):
                data[method].append(float(value))
    return x_values, data


def plot_summary(
    x_values: list[float],
    data: dict[str, list[float]],
    metric: str,
    n: str,
    alpha: str,
    rho: str,
    seed_label: str,
    output_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(6, 4.2))
    plot_methods = COVERAGE_METHODS if metric == "CI_coverage" else METHODS
    for method in plot_methods:
        if method not in data:
            continue
        ax.plot(
            x_values,
            data[method],
            label=method.replace("SAIC", "S-AIC").replace("SBIC", "S-BIC"),
            color=COLORS[method],
            marker=MARKERS[method],
            linewidth=1.8,
            markersize=5,
        )

    if metric == "CI_coverage":
        coverage_values = [
            value
            for method in plot_methods
            if method in data
            for value in data[method]
        ]
        y_min = min(coverage_values)
        ax.axhline(0.95, color="0.45", linestyle="--", linewidth=1)
        ax.set_ylabel("Coverage")
        ax.set_ylim(y_min, 1.0)
    elif metric == "CI_len":
        ax.set_ylabel("Interval length")
    else:
        ax.set_ylabel("Risk")

    ax.set_xlabel(r"$R^2_{\mathrm{signal}}$")
    ax.set_title(f"n={n}, alpha={alpha}, rho={rho}, {seed_label}")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def mean_records(records: list[tuple[list[float], dict[str, list[float]]]]):
    x_values = records[0][0]
    means: dict[str, list[float]] = {}
    for method in records[0][1]:
        columns = [record[1][method] for record in records if method in record[1]]
        means[method] = [sum(values) / len(values) for values in zip(*columns)]
    return x_values, means


def main() -> int:
    input_dir = Path(sys.argv[1]) if len(sys.argv) >= 2 else Path("output/design2")
    plot_dir = Path(sys.argv[2]) if len(sys.argv) >= 3 else input_dir / "plots"
    metric_filter = sys.argv[3] if len(sys.argv) >= 4 else None
    by_seed_dir = plot_dir / "by_seed"
    mean_dir = plot_dir / "seed_average"
    by_seed_dir.mkdir(parents=True, exist_ok=True)
    mean_dir.mkdir(parents=True, exist_ok=True)

    grouped: dict[tuple[str, str, str, str, str], list[tuple[str, list[float], dict[str, list[float]]]]] = defaultdict(list)
    plotted = 0

    for path in sorted(input_dir.glob("*.csv")):
        match = FILE_RE.match(path.name)
        if not match:
            continue
        metric, n, alpha, rho, num_rep, seed = match.groups()
        if metric_filter is not None and metric != metric_filter:
            continue
        x_values, data = read_summary(path)
        output_path = by_seed_dir / f"{metric}_{n}_{alpha}_{rho}_{num_rep}_{seed}.pdf"
        plot_summary(x_values, data, metric, n, alpha, rho, f"seed={seed}", output_path)
        grouped[(metric, n, alpha, rho, num_rep)].append((seed, x_values, data))
        plotted += 1

    for (metric, n, alpha, rho, num_rep), records_with_seed in grouped.items():
        records_with_seed.sort(key=lambda item: int(item[0]))
        seeds = [item[0] for item in records_with_seed]
        records = [(item[1], item[2]) for item in records_with_seed]
        x_values, data = mean_records(records)
        output_path = mean_dir / f"{metric}_{n}_{alpha}_{rho}_{num_rep}_seeds_{seeds[0]}-{seeds[-1]}_mean.pdf"
        plot_summary(x_values, data, metric, n, alpha, rho, f"seeds {seeds[0]}-{seeds[-1]} mean", output_path)
        plotted += 1

    print(f"Plotted {plotted} figures under {plot_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
