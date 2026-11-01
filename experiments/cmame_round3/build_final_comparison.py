"""Build the manuscript comparison from verified Round-2 reruns only."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent/"cmame_round2"/"00_driver_recovery"/"raw_data"/"driver_recovery.json"
OUT = ROOT/"02_final_comparison"


def checkerboard_key(row: dict) -> tuple[str, int]:
    names = {
        "simply_supported_plate": "plate",
        "Scordelis_Lo_roof": "roof",
        "MacNeal_Harder_M3": "M3",
    }
    return names[row["benchmark"]], int(row["n"])


def midedge_key(row: dict) -> tuple[str, int]:
    if row["method"] == "MidedgeTan_linear":
        return "M3", int(row["mesh_cells"])
    if row.get("benchmark") == "Scordelis_Lo_roof":
        return "roof", int(row["mesh_cells_each_direction"])
    return "plate", int(row["mesh_cells"])


def midedge_values(row: dict, benchmark: str, n: int) -> tuple[float, int, int]:
    if benchmark == "plate":
        return (float(row["center_deflection"]), int(row["total_scalar_dof"]),
                int(row["active_linear_dof"]))
    if benchmark == "roof":
        return (float(row["downward_free_edge_midspan_displacement"]),
                int(row["total_scalar_dof"]), int(row["active_dofs"]))
    # The historical M3 wrapper reports its reduced dimension but not the
    # already-known full state formula: 3 vertex states plus one state per
    # triangulation edge = 6*n^2 + 8*n + 3.
    return (float(row["minus_u_B"]), 6*n*n+8*n+3, int(row["active_dof"]))


def build_rows() -> list[dict]:
    payload = json.loads(SOURCE.read_text())
    checkerboard = {checkerboard_key(row): row for row in payload["checkerboard_rows"]}
    midedge = {midedge_key(row): row for row in payload["midedgetan_rows"]}
    references = {"plate": 4.4360891e-4, "roof": 0.3024, "M3": 0.09352155}
    order = {"plate": (8, 16, 32, 64), "roof": (12, 16, 24, 32),
             "M3": (6, 8, 12, 16)}
    rows = []
    for benchmark, meshes in order.items():
        for n in meshes:
            cb = checkerboard[(benchmark, n)]
            mid = midedge[(benchmark, n)]
            mid_response, mid_raw, mid_active = midedge_values(mid, benchmark, n)
            reference = references[benchmark]
            rows.append({
                "benchmark": benchmark, "mesh_n": n, "reference": reference,
                "checkerboard_response": cb["displacement"],
                "checkerboard_error_percent": 100*cb["relative_error"],
                "midedgetan_response": mid_response,
                "midedgetan_error_percent": 100*abs(mid_response-reference)/reference,
                "checkerboard_raw_scalar_states": cb["raw_scalar_states"],
                "midedgetan_raw_scalar_states": mid_raw,
                "checkerboard_active_states": cb["active_states"],
                "midedgetan_active_states": mid_active,
                "raw_state_ratio_checkerboard_over_midedgetan": (
                    cb["raw_scalar_states"]/mid_raw),
                "provenance": "Round-2 recovered manuscript drivers / official LibShell",
            })
    return rows


def write_tables(rows: list[dict]) -> None:
    raw = OUT/"raw_data"; raw.mkdir(parents=True, exist_ok=True)
    with (raw/"final_comparative_evidence.csv").open(
            "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader(); writer.writerows(rows)
    (raw/"final_comparative_evidence.json").write_text(
        json.dumps(rows, indent=2)+"\n", encoding="utf-8")

    header = (
        "| Benchmark | n | Checkerboard response (error) | MidedgeTan response "
        "(error) | Checkerboard raw / active | MidedgeTan raw / active |\n"
        "|---|---:|---:|---:|---:|---:|\n")
    lines = [header]
    for row in rows:
        response_format = ".8g"
        lines.append(
            f"| {row['benchmark']} | {row['mesh_n']} | "
            f"{row['checkerboard_response']:{response_format}} "
            f"({row['checkerboard_error_percent']:.4g}%) | "
            f"{row['midedgetan_response']:{response_format}} "
            f"({row['midedgetan_error_percent']:.4g}%) | "
            f"{row['checkerboard_raw_scalar_states']} / "
            f"{row['checkerboard_active_states']} | "
            f"{row['midedgetan_raw_scalar_states']} / "
            f"{row['midedgetan_active_states']} |\n")
    (OUT/"MANUSCRIPT_COMPARISON_TABLE.md").write_text(
        "# Final comparative evidence table\n\n"+"".join(lines)+
        "\nCPU and RSS are intentionally excluded because the verified methods use "
        "different backends. MITC4 is not included because no provenance-closed "
        "implementation is available.\n", encoding="utf-8")

    latex = [
        "\\begin{tabular}{llrrrrrr}\n",
        "\\toprule\n",
        "Case & $n$ & CB response & CB err. (\\%) & MET response & MET err. "
        "(\\%) & states CB/MET \\\\\n",
        "\\midrule\n",
    ]
    for row in rows:
        latex.append(
            f"{row['benchmark']} & {row['mesh_n']} & "
            f"{row['checkerboard_response']:.7g} & "
            f"{row['checkerboard_error_percent']:.4g} & "
            f"{row['midedgetan_response']:.7g} & "
            f"{row['midedgetan_error_percent']:.4g} & "
            f"{row['checkerboard_raw_scalar_states']}/"
            f"{row['midedgetan_raw_scalar_states']} \\\\\n")
    latex += ["\\bottomrule\n", "\\end{tabular}\n"]
    (OUT/"MANUSCRIPT_COMPARISON_TABLE.tex").write_text("".join(latex), encoding="utf-8")


def make_figure(rows: list[dict]) -> None:
    figures = OUT/"figures"; figures.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 8, "axes.labelsize": 8,
                         "legend.fontsize": 7, "lines.linewidth": 1.2,
                         "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.45))
    for ax, benchmark in zip(axes, ("plate", "roof", "M3")):
        group = [r for r in rows if r["benchmark"] == benchmark]
        ax.loglog([r["checkerboard_raw_scalar_states"] for r in group],
                  [r["checkerboard_error_percent"] for r in group],
                  "o-", color="0.15", label="Checkerboard")
        ax.loglog([r["midedgetan_raw_scalar_states"] for r in group],
                  [r["midedgetan_error_percent"] for r in group],
                  "s--", color="0.60", label="MidedgeTan")
        ax.set_title(benchmark)
        ax.set_xlabel("raw scalar states")
        ax.grid(True, which="both", color="0.88", linewidth=0.5)
    axes[0].set_ylabel("continuum-reference error [%]")
    axes[-1].legend(frameon=False)
    fig.tight_layout()
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(figures/f"Figure_Comparison_accuracy_per_state.{suffix}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)


def run() -> dict:
    rows = build_rows(); write_tables(rows); make_figure(rows)
    summary = {
        "rows": len(rows),
        "maximum_raw_state_ratio": max(
            r["raw_state_ratio_checkerboard_over_midedgetan"] for r in rows),
        "minimum_raw_state_ratio": min(
            r["raw_state_ratio_checkerboard_over_midedgetan"] for r in rows),
        "cpu_or_memory_claim_included": False,
        "mitc4_included": False,
    }
    (OUT/"raw_data"/"summary.json").write_text(
        json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
