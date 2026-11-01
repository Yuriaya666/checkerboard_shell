"""Figures and summary for the aircraft/automotive application studies."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE/"raw_data"
FIGURES = HERE/"figures"
FIGURES.mkdir(parents=True, exist_ok=True)


def load(name):
    return json.loads((DATA/name).read_text())


def aircraft_figure(data):
    xy = np.asarray(data["observation_xy"])
    target = np.asarray(data["target"])
    optimized = np.asarray(data["checkerboard"]["optimized"])
    figure, axes = plt.subplots(1, 2, figsize=(10.5, 3.8))
    for axis, values, title in zip(
        axes, (target, optimized), ("target $u_z$", "optimized $u_z$")):
        scatter = axis.scatter(xy[:, 0], xy[:, 1], c=values,
                               cmap="viridis", s=45)
        figure.colorbar(scatter, ax=axis, label="$u_z$")
        axis.set_xlabel("$x$"); axis.set_ylabel("$y$")
        axis.set_title(title)
    figure.suptitle("Aircraft morphing-skin panel: target vs optimized")
    figure.tight_layout()
    figure.savefig(FIGURES/"Figure_R6_aircraft_skin.png", dpi=220)
    figure.savefig(FIGURES/"Figure_R6_aircraft_skin.pdf")
    plt.close(figure)


def automotive_figure(data):
    target = np.asarray(data["target"])
    optimized = np.asarray(data["checkerboard"]["optimized"])
    index = np.arange(target.size)
    figure, axis = plt.subplots(figsize=(6.8, 3.6))
    axis.plot(index, target, "o-", color="0.55", label="target gap")
    axis.plot(index, optimized, "s-", color="0.1", label="optimized gap")
    axis.set_xlabel("gap measurement point")
    axis.set_ylabel("$u_z$")
    axis.set_title("Automotive panel: gap/flush target")
    axis.grid(True, alpha=0.25)
    axis.legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(FIGURES/"Figure_R6_automotive_panel.png", dpi=220)
    figure.savefig(FIGURES/"Figure_R6_automotive_panel.pdf")
    plt.close(figure)


def cost_figure(cases):
    labels = ["Aircraft skin", "Automotive panel"]
    cases = cases[:2]
    cb_time = [case["checkerboard"]["total_seconds"] for case in cases]
    mt_time = [case["mitc4"]["total_seconds"] for case in cases]
    cb_work = [case["checkerboard"]["solves"] for case in cases]
    mt_work = [case["mitc4"]["analyses"] for case in cases]
    cb_dof = [case["free_dofs"] for case in cases]
    mt_dof = [case["mitc4"].get("free_dofs", 0) for case in cases]
    x = np.arange(len(cases))
    width = 0.35
    figure, axes = plt.subplots(1, 3, figsize=(13.0, 3.8))
    for axis, a, b, title, scale in (
        (axes[0], cb_time, mt_time, "total time [s]", "log"),
        (axes[1], cb_work, mt_work, "work units", "log"),
        (axes[2], cb_dof, mt_dof, "free DOF", "linear"),
    ):
        axis.bar(x - width/2, a, width, color="0.15", label="Checkerboard")
        axis.bar(x + width/2, b, width, color="0.65", label="MITC4")
        axis.set_xticks(x); axis.set_xticklabels(labels, rotation=10)
        axis.set_title(title)
        if scale == "log":
            axis.set_yscale("log")
            axis.set_ylim(bottom=min(a + b) * 0.7, top=max(a + b) * 1.7)
        else:
            axis.set_ylim(0.0, max(a + b) * 1.18)
        axis.legend(fontsize=7)
        for index, value in enumerate(a):
            axis.text(index - width/2, value, f"{value:.2f}",
                      ha="center", va="bottom", fontsize=7)
        for index, value in enumerate(b):
            axis.text(index + width/2, value, f"{value:.2f}",
                      ha="center", va="bottom", fontsize=7)
    figure.tight_layout()
    figure.savefig(FIGURES/"Figure_R6_application_cost.png", dpi=220)
    figure.savefig(FIGURES/"Figure_R6_application_cost.pdf")
    plt.close(figure)


def main():
    aircraft = load("aircraft_skin.json")
    automotive = load("automotive_panel.json")
    aircraft_figure(aircraft)
    automotive_figure(automotive)
    cost_figure([aircraft, automotive])
    print("figures written")


if __name__ == "__main__":
    main()
