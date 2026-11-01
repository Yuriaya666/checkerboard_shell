"""Full 3V plate: convergence figure and cross-check table."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from experiments.cmame_round4.common import DATA, FIGURES, OUT, write_csv


def load_rows() -> list[dict]:
    return json.loads((DATA/"unified_comparison.json").read_text())


def make() -> dict:
    rows = [row for row in load_rows() if row["benchmark"] == "plate"]
    figure, axis = plt.subplots(figsize=(5.4, 4.0))
    for method, marker, style, label, colour in (
        ("Checkerboard-full-3V", "o", "-", "Checkerboard (full 3V)", "0.1"),
        ("Checkerboard-scalar", "^", ":", "Checkerboard (scalar branch)", "0.45"),
        ("MidedgeTan", "s", "--", "MidedgeTan", "0.65"),
    ):
        selected = [row for row in rows if row["method"] == method]
        axis.plot([row["n"] for row in selected],
                  [100.0*float(row["relative_error"]) for row in selected],
                  marker=marker, linestyle=style, color=colour, label=label,
                  markersize=5)
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("cells per direction $n$")
    axis.set_ylabel("center-deflection error [%]")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(fontsize=8)
    figure.tight_layout()
    FIGURES.mkdir(parents=True, exist_ok=True)
    figure.savefig(FIGURES/"Figure_R4_full_plate.png", dpi=220)
    figure.savefig(FIGURES/"Figure_R4_full_plate.pdf")
    plt.close(figure)
    write_csv(DATA/"plate_full_crosscheck.csv", rows)
    return {
        "rows": len(rows),
        "figure": str(FIGURES/"Figure_R4_full_plate.png"),
    }


if __name__ == "__main__":
    print(json.dumps(make(), indent=2))
