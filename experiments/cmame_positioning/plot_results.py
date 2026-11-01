"""Plot the fixed-input positioning diagnostic from saved numerical data."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
FIGURES = HERE / "figures"


def main():
    FIGURES.mkdir(exist_ok=True)
    with (RESULTS / "summary.csv").open(newline="") as file:
        rows = list(csv.DictReader(file))
    styles = {
        "CB-Full": ("#185a9d", "o"),
        "ShellMITC4(Pq)": ("#c15b3e", "s"),
        "ASDShellQ4(Pq)": ("#527d4b", "^"),
    }
    plt.rcParams.update({"font.size": 9, "pdf.fonttype": 42,
                         "svg.fonttype": "none", "axes.spines.top": False,
                         "axes.spines.right": False})

    fig, ax = plt.subplots(figsize=(5.4, 3.7), constrained_layout=True)
    for method, (color, marker) in styles.items():
        subset = sorted((r for r in rows if r["method"] == method and
                         r["response_at_zero_m"]), key=lambda r: int(r["n"]))
        ax.plot([int(r["n"]) for r in subset],
                [-1e3*float(r["response_at_zero_m"]) for r in subset],
                marker=marker, color=color, linewidth=1.4, markersize=4,
                label=method)
    ax.set(xlabel="Cells per side, n", ylabel="Mean loaded-edge deflection (mm)")
    ax.set_xlim(6, 83)
    ax.legend(frameon=False, fontsize=8)
    for extension in ("svg", "pdf"):
        fig.savefig(FIGURES / f"fixed_input_response.{extension}")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.4, 3.7), constrained_layout=True)
    for method, (color, marker) in styles.items():
        subset = sorted((r for r in rows if r["method"] == method and
                         r["central_design_sensitivity_m_per_unit"]),
                        key=lambda r: int(r["n"]))
        ax.plot([int(r["n"]) for r in subset],
                [1e3*float(r["central_design_sensitivity_m_per_unit"])
                 for r in subset], marker=marker, color=color,
                linewidth=1.4, markersize=4, label=method)
    ax.set(xlabel="Cells per side, n",
           ylabel="Deflection sensitivity to geometry parameter (mm/unit)")
    ax.set_xlim(6, 83)
    ax.legend(frameon=False, fontsize=8)
    for extension in ("svg", "pdf"):
        fig.savefig(FIGURES / f"design_sensitivity.{extension}")
    plt.close(fig)


if __name__ == "__main__":
    main()
