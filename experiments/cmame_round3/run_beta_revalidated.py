"""Beta_M sensitivity using only the recovered manuscript benchmark drivers."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg
import torch

from checkerboard_shell.geometry.mesh import CheckerboardMesh
from checkerboard_shell.material.isotropic import IsotropicMaterial
from checkerboard_shell.mechanics.energy import make_targets
from checkerboard_shell.mechanics.low_rank import (
    LowRankTangent, assemble_paper_low_rank_tangent,
)
from experiments.cmame_round2.recovered_drivers import (
    M3_REFERENCE, PLATE_REFERENCE, ROOF_REFERENCE, m3_problem,
    plate_stiffness, roof_problem,
)


ROOT = Path(__file__).resolve().parent
OUT = ROOT/"01_beta_revalidated"
MULTIPLIERS = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0)


def positive_spectrum(matrix: np.ndarray) -> dict:
    eigenvalues = scipy.linalg.eigvalsh(
        0.5*(matrix+matrix.T), check_finite=False, driver="evr")
    scale = max(float(np.max(np.abs(eigenvalues))), 1e-300)
    cutoff = 1e-10*scale
    positive = eigenvalues[eigenvalues > cutoff]
    return {
        "eigenvalue_cutoff": cutoff,
        "smallest_positive_eigenvalues": positive[:5].tolist(),
        "smallest_positive_eigenvalue": (
            float(positive[0]) if positive.size else None),
        "null_or_near_null_count": int(np.count_nonzero(np.abs(eigenvalues) <= cutoff)),
        "negative_count": int(np.count_nonzero(eigenvalues < -cutoff)),
        "maximum_eigenvalue": float(eigenvalues[-1]),
        "condition_estimate_positive": (
            float(positive[-1]/positive[0]) if positive.size else None),
    }


def plate_rows() -> list[dict]:
    n = 16
    stiffness = plate_stiffness(n, include_xw=True)
    grid = np.arange((n+1)**2).reshape(n+1, n+1)
    active = grid[1:-1, 1:-1].reshape(-1)
    reduced = stiffness[active][:, active].toarray()
    force = np.full(active.size, 1.0e-8/(n*n))
    solution = np.linalg.solve(reduced, force)
    observable = float(solution[np.flatnonzero(active == grid[n//2, n//2])[0]])
    spectrum = positive_spectrum(reduced)
    energy = 0.5*float(solution@force)
    rows = []
    for multiplier in MULTIPLIERS:
        rows.append({
            "benchmark": "plate", "n": n, "c_beta": multiplier,
            "status": "SUCCESS", "observable": observable,
            "reference": PLATE_REFERENCE,
            "relative_error": abs(observable-PLATE_REFERENCE)/PLATE_REFERENCE,
            "total_internal_energy_quadratic": energy,
            "xm_energy_quadratic": 0.0, "xm_energy_fraction": 0.0,
            "dxm_channel_status": (
                "NOT_ACTIVE_IN_RECOVERED_SCALAR_TRANSVERSE_BRANCH"),
            "dxm_completion_stiffness": None,
            "dxm_total_directional_stiffness": None,
            "dxm_response_amplitude": None,
            "dxm_gradient_norm": 0.0,
            "blind_direction_interpretation": (
                "not diagnosable in the manuscript scalar transverse branch"),
            **spectrum,
        })
    return rows


def curved_rows(name: str, problem: dict) -> list[dict]:
    reference = problem["reference"]
    tangent, residual, audit = assemble_paper_low_rank_tangent(
        torch.zeros_like(reference), reference, problem["hx"], problem["hy"],
        problem["targets"], problem["material"], variant="CB-Full")
    alpha = tangent.coefficient
    if "active" in problem:
        active = problem["active"]
        local = tangent.local[active][:, active].tocsr()
        g = tangent.vector[active]
        force = problem["force"]
        projector = np.eye(active.size)
        full_from_reduced = None
    else:
        local = tangent.local
        g = tangent.vector
        force = problem["force"]
        constraint = problem["constraints"].toarray()
        projector = scipy.linalg.null_space(constraint, rcond=1e-12)
        full_from_reduced = projector

    local_projected = projector.T @ (local @ projector)
    g_projected = projector.T @ g
    force_projected = projector.T @ force
    gnorm = float(np.linalg.norm(g_projected))
    z = g_projected/gnorm
    rows = []
    for multiplier in MULTIPLIERS:
        coefficient = multiplier*alpha
        matrix = local_projected + coefficient*np.outer(g_projected, g_projected)
        reduced_solution = np.linalg.solve(matrix, force_projected)
        if full_from_reduced is None:
            full = np.zeros(reference.numel())
            full[active] = reduced_solution
        else:
            full = full_from_reduced @ reduced_solution
        if name == "roof":
            observable = -float(full[problem["output_raw"]])
            target = ROOF_REFERENCE
        else:
            observable = -float(problem["point_b"][1]@full)
            target = M3_REFERENCE
        total_energy = 0.5*float(reduced_solution@matrix@reduced_solution)
        xm_energy = 0.5*coefficient*float(g_projected@reduced_solution)**2
        spectrum = positive_spectrum(matrix)
        rows.append({
            "benchmark": name, "n": 16 if name == "roof" else 12,
            "c_beta": multiplier, "status": "SUCCESS",
            "observable": observable, "reference": target,
            "relative_error": abs(observable-target)/target,
            "total_internal_energy_quadratic": total_energy,
            "xm_energy_quadratic": xm_energy,
            "xm_energy_fraction": xm_energy/total_energy if total_energy else 0.0,
            "dxm_channel_status": "ACTIVE",
            "dxm_completion_stiffness": coefficient*gnorm*gnorm,
            "dxm_total_directional_stiffness": float(z@matrix@z),
            "dxm_response_amplitude": float(g_projected@reduced_solution),
            "dxm_gradient_norm": gnorm,
            "reference_residual_inf": float(np.max(np.abs(residual))),
            "rank_one_coefficient_at_c_beta_1": alpha,
            "local_tangent_nnz": int(local.nnz),
            **spectrum,
        })
    return rows


def add_relative_changes(rows: list[dict]) -> None:
    for benchmark in ("plate", "roof", "M3"):
        group = [row for row in rows if row["benchmark"] == benchmark]
        baseline = next(row for row in group if row["c_beta"] == 1.0)["observable"]
        for row in group:
            row["relative_change_vs_c_beta_1"] = (
                abs(row["observable"]-baseline)/abs(baseline))


def write_rows(rows: list[dict]) -> None:
    raw = OUT/"raw_data"; raw.mkdir(parents=True, exist_ok=True)
    keys: list[str] = []
    flat_rows = []
    for row in rows:
        flattened = dict(row)
        flattened["smallest_positive_eigenvalues"] = json.dumps(
            row["smallest_positive_eigenvalues"])
        flat_rows.append(flattened)
        keys.extend(key for key in flattened if key not in keys)
    with (raw/"beta_revalidated.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader(); writer.writerows(flat_rows)
    (raw/"beta_revalidated.json").write_text(
        json.dumps(rows, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def make_figure(rows: list[dict]) -> None:
    figures = OUT/"figures"; figures.mkdir(parents=True, exist_ok=True)
    styles = {
        "plate": dict(color="0.15", marker="o", linestyle="-"),
        "roof": dict(color="0.45", marker="s", linestyle="--"),
        "M3": dict(color="0.70", marker="^", linestyle="-."),
    }
    plt.rcParams.update({"font.size": 8, "axes.labelsize": 8,
                         "legend.fontsize": 7, "lines.linewidth": 1.2,
                         "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.45))
    for benchmark in styles:
        group = [r for r in rows if r["benchmark"] == benchmark]
        axes[0].plot([r["c_beta"] for r in group],
                     [r["relative_change_vs_c_beta_1"] for r in group],
                     label=benchmark, **styles[benchmark])
        axes[1].plot([r["c_beta"] for r in group],
                     [r["xm_energy_fraction"] for r in group],
                     label=benchmark, **styles[benchmark])
    for benchmark in ("roof", "M3"):
        group = [r for r in rows if r["benchmark"] == benchmark]
        axes[2].plot([r["c_beta"] for r in group],
                     [r["dxm_completion_stiffness"] for r in group],
                     label=benchmark, **styles[benchmark])
    axes[0].set_ylabel("relative response change vs $c_\\beta=1$")
    axes[1].set_ylabel("$E_{X_M}/E_{int}$")
    axes[2].set_ylabel("completion-added $DX_M$ stiffness")
    for ax in axes:
        ax.set_xlabel("$c_\\beta=\\beta_M/(8\\mu)$")
        ax.set_xscale("symlog", linthresh=0.25)
        ax.grid(True, which="both", color="0.88", linewidth=0.5)
    axes[0].set_yscale("symlog", linthresh=1e-12)
    axes[2].set_yscale("symlog", linthresh=1.0)
    axes[1].legend(frameon=False)
    axes[2].legend(frameon=False)
    fig.tight_layout()
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(figures/f"Figure_Beta_revalidated.{suffix}", dpi=300,
                    bbox_inches="tight")
    plt.close(fig)


def summarize(rows: list[dict]) -> dict:
    summary = {}
    manuscript_values = {
        "plate": 0.0004502102359571575,
        "roof": 0.3232229477308407,
        "M3": 0.09146271702523975,
    }
    for benchmark in ("plate", "roof", "M3"):
        group = [r for r in rows if r["benchmark"] == benchmark]
        plateau = [r for r in group if 0.25 <= r["c_beta"] <= 4.0]
        c0 = next(r for r in group if r["c_beta"] == 0.0)
        c1 = next(r for r in group if r["c_beta"] == 1.0)
        summary[benchmark] = {
            "maximum_relative_response_change_c_beta_0p25_to_4": max(
                r["relative_change_vs_c_beta_1"] for r in plateau),
            "maximum_relative_response_change_c_beta_0_to_8": max(
                r["relative_change_vs_c_beta_1"] for r in group),
            "xm_energy_fraction_at_c_beta_1": c1["xm_energy_fraction"],
            "null_or_near_null_count_c_beta_0": c0["null_or_near_null_count"],
            "null_or_near_null_count_c_beta_1": c1["null_or_near_null_count"],
            "smallest_positive_eigenvalue_c_beta_0": c0["smallest_positive_eigenvalue"],
            "smallest_positive_eigenvalue_c_beta_1": c1["smallest_positive_eigenvalue"],
            "c_beta_1_manuscript_response_relative_difference": abs(
                c1["observable"]-manuscript_values[benchmark])
                /abs(manuscript_values[benchmark]),
        }
    return summary


def flat_membrane_quotient_diagnostic() -> dict:
    """Isolate the one flat membrane defect without changing the benchmark.

    The recovered plate response is a scalar transverse branch, so beta_M is
    inactive there.  This diagnostic uses the same n=16 flat plate reference
    but restricts the reference tangent to the free in-plane displacement
    block.  Rigid motions and the raw gauge remain in both operators; the
    relevant result is the *change* of nullity when the rank-one completion is
    removed.
    """
    n = 16
    mesh = CheckerboardMesh.rectangle(nx=n, ny=n)
    reference = mesh.reference_vertices
    material = IsotropicMaterial(1.0, 0.3, 0.01)
    targets = make_targets(reference, mesh.hx, mesh.hy)
    tangent, _, _ = assemble_paper_low_rank_tangent(
        torch.zeros_like(reference), reference, mesh.hx, mesh.hy,
        targets, material, variant="CB-Full")
    inplane = np.asarray([3*node+component
                          for node in range((n+1)**2)
                          for component in (0, 1)])
    base = tangent.local[inplane][:, inplane].toarray()
    g = tangent.vector[inplane]
    completed = base+tangent.coefficient*np.outer(g, g)
    values0, vectors0 = scipy.linalg.eigh(
        0.5*(base+base.T), check_finite=False, driver="evr")
    values1 = scipy.linalg.eigvalsh(
        0.5*(completed+completed.T), check_finite=False, driver="evr")
    scale = max(float(np.max(np.abs(values1))), 1e-300)
    cutoff = 1e-10*scale
    null0 = np.abs(values0) <= cutoff
    null1 = np.abs(values1) <= cutoff
    null_basis = vectors0[:, null0]
    overlaps = null_basis.T@g
    target_index = int(np.argmax(np.abs(overlaps)))
    target_mode = null_basis[:, target_index]
    target_added = tangent.coefficient*float(g@target_mode)**2
    return {
        "case": "same recovered plate n=16 reference; free in-plane membrane block",
        "scope": "operator diagnostic, not the scalar transverse benchmark response",
        "eigenvalue_cutoff": cutoff,
        "nullity_c_beta_0": int(np.count_nonzero(null0)),
        "nullity_c_beta_1": int(np.count_nonzero(null1)),
        "additional_nullity_when_completion_removed": int(
            np.count_nonzero(null0)-np.count_nonzero(null1)),
        "target_null_mode_completion_added_stiffness": target_added,
        "target_null_mode_base_stiffness": float(target_mode@base@target_mode),
        "target_null_mode_completed_stiffness": float(
            target_mode@completed@target_mode),
        "smallest_absolute_eigenvalues_c_beta_0": np.sort(np.abs(values0))[:10].tolist(),
        "smallest_absolute_eigenvalues_c_beta_1": np.sort(np.abs(values1))[:10].tolist(),
        "interpretation": (
            "one additional flat in-plane null direction occurs at c_beta=0; "
            "the remaining shared nullity contains rigid and raw-gauge directions"),
    }


def run() -> dict:
    (OUT/"config").mkdir(parents=True, exist_ok=True)
    (OUT/"logs").mkdir(parents=True, exist_ok=True)
    (OUT/"config"/"experiment.json").write_text(json.dumps({
        "c_beta": MULTIPLIERS,
        "definition": "beta_M = c_beta * (8 mu)",
        "drivers": {"plate": 16, "roof": 16, "M3": 12},
        "solution": "recovered manuscript reference-tangent linear benchmark",
        "energy": "quadratic reference-tangent internal energy",
        "eigenvalue_cutoff": "1e-10 times maximum absolute eigenvalue",
        "formulation_changes": "none; c_beta is the sole swept parameter",
    }, indent=2)+"\n", encoding="utf-8")
    rows = plate_rows()
    rows += curved_rows("roof", roof_problem(16))
    rows += curved_rows("M3", m3_problem(12))
    add_relative_changes(rows)
    write_rows(rows); make_figure(rows)
    payload = {"rows": len(rows), "failures": 0,
               "summary": summarize(rows),
               "flat_membrane_quotient_diagnostic": (
                   flat_membrane_quotient_diagnostic())}
    (OUT/"raw_data"/"summary.json").write_text(
        json.dumps(payload, indent=2)+"\n", encoding="utf-8")
    return payload


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
