"""Combine Checkerboard, MidedgeTan, and MITC4 into one comparison."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from experiments.cmame_round4.common import DATA, FIGURES, OUT, write_csv, write_json


def load(name: str) -> list[dict]:
    return json.loads((DATA/name).read_text())


def build_rows() -> list[dict]:
    unified = load("unified_comparison.json")
    mitc4 = load("mitc4.json")
    rows: list[dict] = []

    # Plate: Checkerboard full 3V, MidedgeTan, MITC4.
    for n in (8, 16, 32, 64):
        cb = next(row for row in unified if row["benchmark"] == "plate"
                  and row["n"] == n and row["method"] == "Checkerboard-full-3V")
        mt = next(row for row in unified if row["benchmark"] == "plate"
                  and row["n"] == n and row["method"] == "MidedgeTan")
        mi = next(row for row in mitc4 if row["benchmark"] == "plate"
                  and row["n"] == n)
        rows.append({"benchmark": "plate", "n": n, "method": "Checkerboard",
                     "response": cb["response"], "relative_error": cb["relative_error"],
                     "raw_scalar_states": cb["raw_scalar_states"],
                     "active_states": cb["active_states"]})
        rows.append({"benchmark": "plate", "n": n, "method": "MidedgeTan",
                     "response": mt["response"], "relative_error": mt["relative_error"],
                     "raw_scalar_states": mt["raw_scalar_states"],
                     "active_states": mt["active_states"]})
        rows.append({"benchmark": "plate", "n": n, "method": "MITC4",
                     "response": mi["response"], "relative_error": mi["relative_error"],
                     "raw_scalar_states": mi["raw_scalar_states"],
                     "active_states": mi["active_states"],
                     "notes": mi["notes"]})

    # Roof: historical reference for the main comparison; modern in the table.
    for n in (12, 16, 24, 32):
        cb = next(row for row in unified if row["benchmark"] == "roof"
                  and row["n"] == n and row["method"] == "Checkerboard")
        mt = next(row for row in unified if row["benchmark"] == "roof"
                  and row["n"] == n and row["method"] == "MidedgeTan")
        mi = next(row for row in mitc4 if row["benchmark"] == "roof"
                  and row["n"] == n)
        for method, row in (("Checkerboard", cb), ("MidedgeTan", mt),
                            ("MITC4", mi)):
            rows.append({
                "benchmark": "roof", "n": n, "method": method,
                "response": row["response"],
                "relative_error": row["relative_error"],
                "relative_error_modern": row.get("relative_error_modern"),
                "raw_scalar_states": row["raw_scalar_states"],
                "active_states": row["active_states"],
                "notes": row.get("notes", ""),
            })

    # M3: Checkerboard raw (main) plus midpoint, MidedgeTan raw, MITC4 raw.
    for n in (6, 8, 12, 16):
        cb = next(row for row in unified if row["benchmark"] == "M3"
                  and row["n"] == n and row["method"] == "Checkerboard")
        mt = next(row for row in unified if row["benchmark"] == "M3"
                  and row["n"] == n and row["method"] == "MidedgeTan")
        mi = next(row for row in mitc4 if row["benchmark"] == "M3"
                  and row.get("load_mode") == "raw" and row["n"] == n)
        rows.append({
            "benchmark": "M3", "n": n, "method": "Checkerboard",
            "response": cb["response"], "relative_error": cb["relative_error"],
            "response_midpoint": cb.get("response_midpoint"),
            "relative_error_midpoint": cb.get("relative_error_midpoint"),
            "raw_scalar_states": cb["raw_scalar_states"],
            "active_states": cb["active_states"],
        })
        rows.append({
            "benchmark": "M3", "n": n, "method": "MidedgeTan",
            "response": mt["response"], "relative_error": mt["relative_error"],
            "raw_scalar_states": mt["raw_scalar_states"],
            "active_states": mt["active_states"],
        })
        rows.append({
            "benchmark": "M3", "n": n, "method": "MITC4",
            "response": mi["response"], "relative_error": mi["relative_error"],
            "raw_scalar_states": mi["raw_scalar_states"],
            "active_states": mi["active_states"],
            "notes": mi["notes"],
        })
    # MITC4 M3 convergence points beyond the manuscript mesh range.
    for n in (24, 32, 48, 64):
        mi = next(row for row in mitc4 if row["benchmark"] == "M3"
                  and row.get("load_mode") == "raw" and row["n"] == n)
        rows.append({
            "benchmark": "M3", "n": n, "method": "MITC4",
            "response": mi["response"], "relative_error": mi["relative_error"],
            "raw_scalar_states": mi["raw_scalar_states"],
            "active_states": mi["active_states"],
            "notes": mi["notes"] + "; convergence point beyond manuscript range",
        })
    return rows


def markdown(rows: list[dict]) -> str:
    lines = [
        "# Combined Checkerboard / MidedgeTan / MITC4 comparison",
        "",
        "MITC4: OpenSeesPy `ShellMITC4` with `ElasticMembranePlateSection`, "
        "6 DOF per node.  MidedgeTan: official LibShell "
        "`MidedgeAngleTanFormulation`.  Checkerboard: current production "
        "energy.  Roof errors are historical / modern.",
        "",
        "| benchmark | n | Checkerboard | MidedgeTan | MITC4 | CB / MT / MITC4 raw states |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for benchmark in ("plate", "roof", "M3"):
        meshes = sorted({row["n"] for row in rows if row["benchmark"] == benchmark})
        for n in meshes:
            def get(method):
                return next((row for row in rows if row["benchmark"] == benchmark
                             and row["n"] == n and row["method"] == method), None)
            cb, mt, mi = get("Checkerboard"), get("MidedgeTan"), get("MITC4")
            if mi is None:
                continue
            if cb is None or mt is None:
                lines.append(
                    f'| {benchmark} | {n} | — | — | '
                    f'{100.0*mi["relative_error"]:.3f}% | '
                    f'— / — / {mi["raw_scalar_states"]} |')
                continue
            def text(row):
                if benchmark == "roof":
                    return (f'{100.0*row["relative_error"]:.3f}% / '
                            f'{100.0*row["relative_error_modern"]:.3f}%')
                if benchmark == "M3" and row["method"] == "Checkerboard":
                    return (f'{100.0*row["relative_error"]:.3f}% / '
                            f'{100.0*row["relative_error_midpoint"]:.3f}%')
                return f'{100.0*row["relative_error"]:.3f}%'
            lines.append(
                f'| {benchmark} | {n} | {text(cb)} | {text(mt)} | {text(mi)} | '
                f'{cb["raw_scalar_states"]} / {mt["raw_scalar_states"]} / '
                f'{mi["raw_scalar_states"]} |')
    lines.extend([
        "",
        "M3 Checkerboard column: raw equator vertex / physical-midpoint "
        "reconstruction.  MidedgeTan and MITC4 use the raw equator vertex.  "
        "MITC4 M3 n=24,32,48,64 are convergence points beyond the manuscript range.",
        "",
        "MITC4 is a Mindlin-Reissner shell; at t/a=0.01 its plate response has "
        "a ~0.5% transverse-shear floor relative to the thin-plate reference.  "
        "At t/a=0.001 the error drops to 0.03%/0.001%, confirming the source.",
    ])
    return "\n".join(lines) + "\n"


def figure(rows: list[dict]) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(13.2, 3.9))
    styles = {
        "Checkerboard": ("o", "-", "0.10"),
        "MidedgeTan": ("s", "--", "0.55"),
        "MITC4": ("^", "-.", "0.30"),
    }
    for axis, benchmark in zip(axes, ("plate", "roof", "M3")):
        for method, (marker, style, colour) in styles.items():
            selected = [row for row in rows if row["benchmark"] == benchmark
                        and row["method"] == method]
            if not selected:
                continue
            label = ("MITC4 (Mindlin)" if method == "MITC4"
                     and benchmark == "plate" else method)
            axis.plot([row["raw_scalar_states"] for row in selected],
                      [100.0*float(row["relative_error"]) for row in selected],
                      marker=marker, linestyle=style, color=colour,
                      label=label, markersize=5)
        if benchmark == "M3":
            selected = [row for row in rows if row["benchmark"] == "M3"
                        and row["method"] == "Checkerboard"]
            axis.plot([row["raw_scalar_states"] for row in selected],
                      [100.0*float(row["relative_error_midpoint"]) for row in selected],
                      marker="^", linestyle=":", color="0.45",
                      label="Checkerboard midpoint", markersize=4)
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_title(benchmark)
        axis.set_xlabel("raw scalar states")
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(fontsize=7)
    axes[0].set_ylabel("continuum-reference error [%]")
    figure.tight_layout()
    FIGURES.mkdir(parents=True, exist_ok=True)
    figure.savefig(FIGURES/"Figure_R4_combined_three_methods.png", dpi=220)
    figure.savefig(FIGURES/"Figure_R4_combined_three_methods.pdf")
    plt.close(figure)


def run() -> dict:
    rows = build_rows()
    write_csv(DATA/"combined_three_methods.csv", rows)
    write_json(DATA/"combined_three_methods.json", rows)
    (OUT/"COMBINED_COMPARISON_TABLE.md").write_text(markdown(rows),
                                                    encoding="utf-8")
    figure(rows)
    return {"rows": len(rows)}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
