"""Create publication-size Round-2 figures and summary tables from raw data."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
BOUNDARY = ROOT/"01_sparse_boundary"
SCALING = ROOT/"03_scaling"
DISTORTION = ROOT/"04_distortion_revalidated"


plt.rcParams.update({
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "legend.fontsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "lines.linewidth": 1.2,
    "lines.markersize": 4,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

STYLES = {
    "CB-0": dict(color="0.15", marker="o", linestyle="-"),
    "CB-M": dict(color="0.35", marker="s", linestyle="--"),
    "CB-Full": dict(color="0.55", marker="^", linestyle="-."),
    "MidedgeTan": dict(color="0.05", marker="x", linestyle=":"),
}


def save_all(fig, directory: Path, stem: str) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(directory/f"{stem}.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def scaling_figures(rows: list[dict]) -> None:
    figures = SCALING/"figures"
    cb_methods = ("CB-0", "CB-M", "CB-Full")

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7))
    for ax, benchmark in zip(axes, ("plate", "M3")):
        for method in cb_methods:
            selected = sorted((r for r in rows if r["benchmark"] == benchmark
                               and r["method"] == method), key=lambda r: r["n"])
            states = [r["raw_scalar_states"] for r in selected]
            total = [r["assembly_time_median_s"]
                     + (r.get("factorization_time_median_s") or 0.0)
                     + (r.get("linear_solve_time_median_s") or 0.0)
                     for r in selected]
            ax.loglog(states, total, label=method, **STYLES[method])
        ax.set_title(benchmark)
        ax.set_xlabel("raw scalar states")
        ax.grid(True, which="both", color="0.88", linewidth=0.5)
    axes[0].set_ylabel("assembly + factor/solve time [s]")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    save_all(fig, figures, "Figure_C1_time_scaling")

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7))
    for ax, benchmark in zip(axes, ("plate", "M3")):
        for method in cb_methods + ("MidedgeTan",):
            selected = sorted((r for r in rows if r["benchmark"] == benchmark
                               and r["method"] == method), key=lambda r: r["n"])
            memory = [(r.get("child_peak_rss_mb") if method == "MidedgeTan"
                       else r.get("peak_rss_mb")) for r in selected]
            ax.semilogx([r["raw_scalar_states"] for r in selected],
                        memory,
                        label=method, **STYLES[method])
        ax.set_title(benchmark)
        ax.set_xlabel("raw scalar states")
        ax.grid(True, which="both", color="0.88", linewidth=0.5)
    axes[0].set_ylabel("measured worker peak RSS [MiB]")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    save_all(fig, figures, "Figure_C2_memory_scaling")

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7))
    for ax, benchmark in zip(axes, ("plate", "M3")):
        methods = cb_methods + (("MidedgeTan",) if benchmark == "plate" else ())
        for method in methods:
            selected = sorted((r for r in rows if r["benchmark"] == benchmark
                               and r["method"] == method
                               and r.get("augmented_system_nnz") is not None),
                              key=lambda r: r["n"])
            if selected:
                ax.loglog([r["raw_scalar_states"] for r in selected],
                          [r["augmented_system_nnz"] for r in selected],
                          label=method, **STYLES[method])
        ax.set_title(benchmark)
        ax.set_xlabel("raw scalar states")
        ax.grid(True, which="both", color="0.88", linewidth=0.5)
    axes[0].set_ylabel("augmented-system nnz")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    save_all(fig, figures, "Figure_C3_nnz_scaling")

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7))
    for ax, benchmark in zip(axes, ("plate", "M3")):
        for method in ("CB-Full", "MidedgeTan"):
            selected = sorted((r for r in rows if r["benchmark"] == benchmark
                               and r["method"] == method),
                              key=lambda r: r["raw_scalar_states"])
            ax.loglog([r["raw_scalar_states"] for r in selected],
                      [100*r["relative_error"] for r in selected],
                      label=method, **STYLES[method])
        ax.set_title(benchmark)
        ax.set_xlabel("raw scalar states")
        ax.grid(True, which="both", color="0.88", linewidth=0.5)
    axes[0].set_ylabel("relative observable error [%]")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    save_all(fig, figures, "Figure_D3_error_vs_state_count")


def boundary_figure() -> None:
    legacy = json.loads(
        (BOUNDARY/"raw_data"/"legacy_dense_boundary.json").read_text())
    local = json.loads(
        (BOUNDARY/"raw_data"/"boundary_storage.json").read_text())
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7))
    axes[0].loglog([r["n"] for r in legacy[:4]],
                   [r["old_dense_bytes_measured"]/2**20 for r in legacy[:4]],
                   color="0.15", marker="o", label="old dense (measured)")
    axes[0].loglog([legacy[-1]["n"]],
                   [legacy[-1]["old_dense_bytes_projected"]/2**20],
                   color="0.15", marker="o", fillstyle="none",
                   label="old dense (projected)")
    axes[0].loglog([r["n"] for r in local],
                   [r["new_stencil_bytes_measured"]/2**20 for r in local],
                   color="0.55", marker="s", linestyle="--",
                   label="local stencil (measured)")
    axes[0].set(xlabel="cells per direction n", ylabel="operator storage [MiB]")
    axes[0].grid(True, which="both", color="0.88", linewidth=0.5)
    axes[0].legend(frameon=False)
    axes[1].loglog([r["n"] for r in legacy[:4]],
                   [r["old_geometry_time_median_seconds"] for r in legacy[:4]],
                   color="0.15", marker="o", label="old dense")
    axes[1].loglog([r["n"] for r in local],
                   [r["new_geometry_time_median_seconds"] for r in local],
                   color="0.55", marker="s", linestyle="--", label="local stencil")
    axes[1].set(xlabel="cells per direction n", ylabel="geometry evaluation [s]")
    axes[1].grid(True, which="both", color="0.88", linewidth=0.5)
    axes[1].legend(frameon=False)
    fig.tight_layout()
    save_all(fig, BOUNDARY/"figures", "Figure_S1_boundary_storage")


def plot_grid(ax, values: np.ndarray, **style) -> None:
    for i in range(values.shape[0]):
        ax.plot(*values[i, :, :].T, **style)
    for j in range(values.shape[1]):
        ax.plot(*values[:, j, :].T, **style)


def representative_distorted_meshes() -> None:
    fig = plt.figure(figsize=(7.0, 5.0))
    cases = (("plate", "random"), ("roof", "random"), ("M3", "random"))
    for index, (benchmark, family) in enumerate(cases, 1):
        archive = np.load(DISTORTION/"raw_data"/f"{benchmark}_{family}_030.npz")
        reference = archive["reference"]
        deformed = archive["deformed"]
        displacement = deformed-reference
        characteristic = np.linalg.norm(
            reference.reshape(-1, 3).ptp(axis=0))
        maximum = np.linalg.norm(displacement.reshape(-1, 3), axis=1).max()
        scale = 0.15*characteristic/maximum if maximum > 0 else 1.0
        amplified = reference+scale*displacement
        ax = fig.add_subplot(2, 2, index, projection="3d")
        plot_grid(ax, reference, color="0.72", linewidth=0.45)
        plot_grid(ax, amplified, color="0.08", linewidth=0.65)
        ax.set_title(f"{benchmark}, random δ=0.30 (deformation ×{scale:.1g})")
        ax.set_axis_off()
        ax.set_box_aspect(np.maximum(reference.reshape(-1, 3).ptp(axis=0), 0.2))
    legend = fig.add_subplot(2, 2, 4)
    legend.plot([], [], color="0.72", label="distorted reference grid")
    legend.plot([], [], color="0.08", label="amplified deformed grid")
    legend.legend(loc="center", frameon=False)
    legend.set_axis_off()
    fig.tight_layout()
    save_all(fig, DISTORTION/"figures", "Figure_D2_distorted_meshes")


def comparative_distortion() -> None:
    checkerboard = json.loads(
        (DISTORTION/"raw_data"/"distortion.json").read_text())
    midedge = json.loads(
        (DISTORTION/"raw_data"/"midedgetan_distortion.json").read_text())
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.55), sharey=False)
    kind_style = {"smooth": "-", "random": "--"}
    method_marker = {"CB-Full": "o", "MidedgeTan": "x"}
    for ax, benchmark in zip(axes, ("plate", "roof", "M3")):
        for method, source in (("CB-Full", checkerboard), ("MidedgeTan", midedge)):
            for kind in ("smooth", "random"):
                selected = [r for r in source if r["benchmark"] == benchmark
                            and r["distortion_kind"] == kind
                            and r["status"] == "SUCCESS"]
                ax.plot([r["delta"] for r in selected],
                        [100*r["relative_error"] for r in selected],
                        color=("0.15" if method == "CB-Full" else "0.55"),
                        marker=method_marker[method], linestyle=kind_style[kind],
                        label=f"{method}, {kind}")
        ax.set_title(benchmark)
        ax.set_xlabel("redistribution amplitude δ")
        ax.grid(True, color="0.88", linewidth=0.5)
    axes[0].set_ylabel("relative observable error [%]")
    axes[-1].legend(frameon=False, fontsize=6)
    fig.tight_layout()
    save_all(fig, DISTORTION/"figures", "Figure_D1_distortion_comparison")


def summary_table(rows: list[dict]) -> None:
    selected = []
    for benchmark, n in (("plate", 128), ("M3", 128)):
        for method in ("CB-0", "CB-M", "CB-Full", "MidedgeTan"):
            row = next(r for r in rows if r["benchmark"] == benchmark
                       and r["n"] == n and r["method"] == method)
            selected.append({
                "benchmark": benchmark, "n": n, "method": method,
                "raw_states": row["raw_scalar_states"],
                "active_states": row["active_states"],
                "local_nnz": row.get("local_tangent_nnz"),
                "augmented_nnz": row.get("augmented_system_nnz"),
                "rank_one_storage": row.get("rank_one_storage"),
                "peak_rss_mib": (row.get("child_peak_rss_mb")
                                 if method == "MidedgeTan"
                                 else row.get("peak_rss_mb")),
                "observable": row["observable"],
                "relative_error_percent": 100*row["relative_error"],
            })
    path = SCALING/"raw_data"/"Table_C1_summary.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=selected[0].keys())
        writer.writeheader(); writer.writerows(selected)


def run() -> None:
    rows = json.loads((SCALING/"raw_data"/"scaling.json").read_text())
    boundary_figure()
    scaling_figures(rows)
    comparative_distortion()
    representative_distorted_meshes()
    summary_table(rows)


if __name__ == "__main__":
    run()
