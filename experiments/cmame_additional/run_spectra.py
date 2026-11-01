#!/usr/bin/env python3
"""Run curved-reference and free-boundary spectral experiments."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys
import time
import traceback

import numpy as np
import torch
from scipy.linalg import eigh as dense_eigh
from scipy.sparse.linalg import ArpackNoConvergence, LinearOperator, eigsh

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(HERE))

from common import (  # noqa: E402
    DTYPE, ROOT, SEED, alternating_vector, curved_patch_mesh, ensure_tree,
    formulation_energy, method_unavailable_rows, plot_style, save_figure,
    strip_mesh, write_csv, write_json,
)
from checkerboard_shell import BoundaryConditions, IsotropicMaterial  # noqa: E402
from checkerboard_shell.geometry.core import edge_midpoints  # noqa: E402
from checkerboard_shell.mechanics.energy import make_targets, paper_components  # noqa: E402


torch.set_default_dtype(DTYPE)
torch.manual_seed(SEED)
np.random.seed(SEED)


def log_line(handle, message):
    text = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}"
    print(text, flush=True); handle.write(text+"\n"); handle.flush()


def curved_boundary(mesh):
    mask = torch.zeros_like(mesh.reference_vertices, dtype=torch.bool)
    # Two raw layers give an unambiguous essential clamp and remove both rigid
    # and raw-gauge modes from this *physical blind-mode* experiment.
    mask[:2] = True; mask[-2:] = True; mask[:, :2] = True; mask[:, -2:] = True
    return BoundaryConditions(mask)


def free_strip_boundary(mesh):
    mask = torch.zeros_like(mesh.reference_vertices, dtype=torch.bool)
    # Fixing just one raw node is insufficient: a rigid translation plus an
    # opposite raw gauge can vanish there.  A 2x2 raw corner block removes all
    # rigid/gauge combinations (three nodes still leave one mixed combination).
    # The constrained physical area shrinks under refinement and both target
    # longitudinal edges are otherwise free.
    mask[:2, :2, :] = True
    return BoundaryConditions(mask)


def observation_energy(mesh, material, bc, operator):
    ref = mesh.reference_vertices
    targets = make_targets(ref, mesh.hx, mesh.hy)
    q0 = bc.initial_coordinates(ref)
    boundary = targets.checkerboard.boundary

    def function(q):
        positions = bc.positions(ref, q)
        parts = paper_components(positions, ref, mesh.hx, mesh.hy, targets, material)
        geom = parts["geometry"]
        if operator in ("J_B", "J_B+DXM"):
            change = geom.bulk.first_form-boundary.first_form_b
            value = 0.5*(boundary.area_b[..., None, None]*change.square()).sum()
            if operator == "J_B+DXM": value = value+0.5*parts["X_M"].square()
        elif operator in ("J_bW", "J_bW+DXWrel"):
            change = geom.second_form-boundary.second_form
            weights = boundary.dual_fraction*boundary.area_w
            value = 0.5*(weights[..., None, None]*change.square()).sum()
            if operator == "J_bW+DXWrel":
                value = value+0.5*(boundary.gamma_w[1:-1, 1:-1, None]
                                   * parts["X_W_rel"][1:-1, 1:-1].square()).sum()
        elif operator == "full_tangent":
            value = parts["total"]
        else:
            raise ValueError(operator)
        return value
    return function, q0


def smallest_spectrum(function, q0, k=20):
    n = q0.numel()
    k = min(k, n-2)
    started = time.perf_counter()
    # The production implementation is dense already.  Assembling the exact
    # AD Hessian once is both substantially faster and a more faithful audit
    # than thousands of repeated HVP traces through the dense reconstruction
    # matrices.  ARPACK is then used only as the dense symmetric eigensolver.
    tangent = torch.autograd.functional.hessian(function, q0, vectorize=True)
    tangent = 0.5*(tangent+tangent.T)
    matrix = tangent.detach().cpu().numpy()
    # LAPACK's subset driver is robust for the large exact nullspaces of the
    # observation operators; ARPACK ``which=SM`` can stagnate on those zeros.
    values, vectors = dense_eigh(
        matrix, subset_by_index=[0, k-1], driver="evr", check_finite=False)
    largest = float(dense_eigh(
        matrix, subset_by_index=[n-1, n-1], driver="evr",
        eigvals_only=True, check_finite=False)[0])
    status = "CONVERGED_DENSE_SUBSET"
    return values, vectors, abs(largest), status, 0, time.perf_counter()-started


def mode_metrics(mesh, bc, vector):
    field = torch.zeros_like(mesh.reference_vertices)
    field[bc.free_mask] = torch.as_tensor(vector, dtype=DTYPE)
    mids = edge_midpoints(field)
    physical = math.sqrt(float(torch.linalg.norm(mids.horizontal)**2
                               + torch.linalg.norm(mids.vertical)**2))
    parity = alternating_vector(mesh.reference_vertices.shape)
    alternating = float(torch.linalg.norm((parity*field).sum((0, 1)))
                        / (math.sqrt(mesh.num_vertices)*torch.linalg.norm(field)))
    boundary_mask = torch.zeros(mesh.reference_vertices.shape[:2], dtype=torch.bool)
    boundary_mask[:2] = True; boundary_mask[-2:] = True
    boundary_mask[:, :2] = True; boundary_mask[:, -2:] = True
    boundary_ratio = float(torch.linalg.norm(field[boundary_mask])/torch.linalg.norm(field))
    return field, physical, alternating, boundary_ratio


def classify(value, scale, physical, alternating, boundary_ratio):
    near = abs(value) <= max(scale*1e-9, 1e-14)
    if physical <= 1e-8: return "gauge"
    if near and alternating >= 0.45 and boundary_ratio >= 0.7: return "boundary-alternating"
    if near and alternating >= 0.45: return "alternating/checkerboard-like"
    if near and boundary_ratio >= 0.8: return "boundary-localized"
    if near: return "smooth physical mode"
    return "positive mode"


def analyze_case(mesh, material, bc, operator, context, archive_dir, log):
    function, q0 = observation_energy(mesh, material, bc, operator)
    values, vectors, scale, status, calls, runtime = smallest_spectrum(function, q0)
    rows = []
    fields = []
    for index, value in enumerate(values):
        field, physical, alternating, boundary_ratio = mode_metrics(
            mesh, bc, vectors[:, index])
        rows.append({
            **context, "operator": operator, "mode_index": index,
            "eigenvalue": float(value),
            "singular_value": math.sqrt(max(float(value), 0.0)),
            "largest_abs_eigenvalue": scale,
            "relative_eigenvalue": float(value)/scale if scale else math.nan,
            "physical_midpoint_norm": physical,
            "alternating_score": alternating,
            "boundary_localization_ratio": boundary_ratio,
            "classification": classify(float(value), scale, physical,
                                       alternating, boundary_ratio),
            "solver_status": status, "operator_calls": calls,
            "spectrum_runtime_s": runtime,
            "active_state_count": q0.numel(),
        })
        if index < 6: fields.append(field.numpy())
    archive_dir.mkdir(parents=True, exist_ok=True)
    tag = "_".join(str(context[key]) for key in context)+f"_{operator}"
    np.savez_compressed(archive_dir/f"modes_{tag}.npz",
                        reference=mesh.reference_vertices.numpy(),
                        modes=np.asarray(fields), eigenvalues=values,
                        fixed_mask=bc.fixed_mask.numpy(), quads=mesh.quads)
    log_line(log, f"spectrum {tag}: {status}, calls={calls}, runtime={runtime:.2f}s, "
                  f"lambda0={values[0]:.3e}")
    return rows


def run_curved(log):
    base = ROOT/"05_curved_spectrum"
    sizes = [8, 16, 32]
    operators = ["J_B", "J_B+DXM", "J_bW", "J_bW+DXWrel", "full_tangent"]
    write_json(base/"config"/"experiment.json", {
        "curvature_types": ["sphere", "saddle"], "n": sizes,
        "operators": operators, "requested_modes": 20,
        "BC": "outer two raw layers fixed in all components",
        "near_zero_relative_threshold": 1e-9,
    })
    rows = []
    material = IsotropicMaterial(young=1.0, poisson=0.3, thickness=0.05)
    for kind in ("sphere", "saddle"):
        for n in sizes:
            mesh = curved_patch_mesh(kind, n); bc = curved_boundary(mesh)
            for operator in operators:
                try:
                    rows += analyze_case(mesh, material, bc, operator,
                                         {"curvature": kind, "n": n},
                                         base/"raw_data"/"mode_shapes", log)
                except Exception as error:
                    rows.append({"curvature": kind, "n": n, "operator": operator,
                                 "solver_status": "FAILED", "reason": repr(error)})
                    log_line(log, f"curved {kind} n={n} {operator} failed: {error}")
                    traceback.print_exc(file=log)
    write_csv(base/"raw_data"/"spectra.csv", rows)
    write_json(base/"raw_data"/"spectra.json", rows)
    _plot_curved(rows)
    return rows


def _plot_curved(rows):
    plot_style(); import matplotlib.pyplot as plt
    valid = [r for r in rows if "mode_index" in r]
    styles = {"J_B": "o-", "J_B+DXM": "s--", "J_bW": "^-",
              "J_bW+DXWrel": "d--", "full_tangent": "x-"}
    for kind, figure in (("sphere", "Figure_K1"), ("saddle", "Figure_K2")):
        fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.4), sharey=True)
        for ax, n in zip(axes, [8, 16, 32]):
            for operator, style in styles.items():
                subset = [r for r in valid if r["curvature"] == kind and r["n"] == n
                          and r["operator"] == operator]
                ax.semilogy([r["mode_index"]+1 for r in subset],
                            [max(abs(r["relative_eigenvalue"]), 1e-18) for r in subset],
                            style, ms=2.5, label=operator)
            ax.set_title(f"n={n}"); ax.set_xlabel("mode index"); ax.grid(True, alpha=0.25)
        axes[0].set_ylabel(r"$|\lambda|/\lambda_{max}$")
        axes[-1].legend(frameon=False, fontsize=8)
        save_figure(fig, ROOT/"05_curved_spectrum"/"figures"/figure); plt.close(fig)
    # Classification diagnostic doubles as representative soft-mode summary.
    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    for kind, marker in (("sphere", "o"), ("saddle", "s")):
        subset = [r for r in valid if r["curvature"] == kind and r["operator"] == "full_tangent"]
        ax.scatter([r["alternating_score"] for r in subset],
                   [max(abs(r["relative_eigenvalue"]), 1e-18) for r in subset],
                   marker=marker, s=12, facecolors="none", label=kind)
    ax.set(xlabel="normalized alternating score", ylabel=r"$|\lambda|/\lambda_{max}$", yscale="log")
    ax.grid(True, alpha=0.25); ax.legend(frameon=False)
    save_figure(fig, ROOT/"05_curved_spectrum"/"figures"/"Figure_K3"); plt.close(fig)


def run_boundary(log):
    base = ROOT/"06_boundary_modes"
    sizes = [8, 16, 32]
    operators = ["J_bW", "J_bW+DXWrel", "full_tangent"]
    write_json(base/"config"/"experiment.json", {
        "patches": ["flat", "cylinder"], "n": sizes, "operators": operators,
        "free_edges": "both longitudinal y-boundaries remain unconstrained",
        "rigid_and_gauge_removal": "2x2 corner raw block; removes all mixed rigid/gauge combinations",
        "boundary_ratio": "two-layer displacement norm / total norm",
    })
    rows = []
    material = IsotropicMaterial(young=1.0, poisson=0.3, thickness=0.05)
    for kind in ("flat", "cylinder"):
        for n in sizes:
            mesh = strip_mesh(kind, n); bc = free_strip_boundary(mesh)
            for operator in operators:
                try:
                    rows += analyze_case(mesh, material, bc, operator,
                                         {"patch": kind, "n": n},
                                         base/"raw_data"/"mode_shapes", log)
                except Exception as error:
                    rows.append({"patch": kind, "n": n, "operator": operator,
                                 "solver_status": "FAILED", "reason": repr(error)})
                    log_line(log, f"boundary {kind} n={n} {operator} failed: {error}")
                    traceback.print_exc(file=log)
    contexts = [{"patch": kind, "n": n, "operator": "low_spectrum"}
                for kind in ("flat", "cylinder") for n in sizes]
    rows += method_unavailable_rows(
        "MidedgeTan", contexts,
        "MidedgeTan is absent from the repository; no unvalidated spectrum substituted.")
    rows += method_unavailable_rows(
        "MITC4", contexts,
        "A verified MITC4 backend is absent; no unvalidated spectrum substituted.")
    write_csv(base/"raw_data"/"spectra.csv", rows)
    write_json(base/"raw_data"/"spectra.json", rows)
    _plot_boundary(rows)
    return rows


def _plot_boundary(rows):
    plot_style(); import matplotlib.pyplot as plt
    valid = [r for r in rows if "mode_index" in r]
    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    for kind, marker in (("flat", "o"), ("cylinder", "s")):
        x, y = [], []
        for n in [8, 16, 32]:
            subset = [r for r in valid if r["patch"] == kind and r["n"] == n
                      and r["operator"] == "full_tangent"]
            if subset:
                positive = [abs(r["eigenvalue"]) for r in subset if r["eigenvalue"] > 0]
                if positive: x.append(n); y.append(min(positive))
        ax.loglog(x, y, marker=marker, label=kind)
    ax.set(xlabel="n", ylabel="smallest reported positive eigenvalue")
    ax.grid(True, which="both", alpha=0.25); ax.legend(frameon=False)
    save_figure(fig, ROOT/"06_boundary_modes"/"figures"/"Figure_E1"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    for kind, marker in (("flat", "o"), ("cylinder", "s")):
        subset = [r for r in valid if r["patch"] == kind and r["operator"] == "full_tangent"]
        ax.scatter([max(abs(r["relative_eigenvalue"]), 1e-18) for r in subset],
                   [r["boundary_localization_ratio"] for r in subset],
                   marker=marker, s=12, facecolors="none", label=kind)
    ax.set(xlabel=r"$|\lambda|/\lambda_{max}$", ylabel="boundary localization", xscale="log")
    ax.grid(True, alpha=0.25); ax.legend(frameon=False)
    save_figure(fig, ROOT/"06_boundary_modes"/"figures"/"Figure_E2"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    subset = [r for r in valid if r["operator"] == "full_tangent"]
    ax.scatter([r["alternating_score"] for r in subset],
               [r["boundary_localization_ratio"] for r in subset],
               c=[math.log10(max(abs(r["relative_eigenvalue"]), 1e-18)) for r in subset],
               cmap="Greys", s=14)
    ax.set(xlabel="alternating score", ylabel="boundary localization")
    ax.grid(True, alpha=0.25)
    save_figure(fig, ROOT/"06_boundary_modes"/"figures"/"Figure_E3"); plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--only", choices=["curved", "boundary"])
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    torch.set_num_threads(args.threads); ensure_tree()
    with (ROOT/"spectral_run.log").open("a") as log:
        log_line(log, f"Spectral run start, only={args.only}, threads={args.threads}")
        if args.only in (None, "curved"): run_curved(log)
        if args.only in (None, "boundary"): run_boundary(log)
        log_line(log, "Spectral run finished")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
