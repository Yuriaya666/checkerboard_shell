#!/usr/bin/env python3
"""Run CMAME experiments 1--4 and cross-cutting consistency checks."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
import time
import traceback

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(HERE))

from common import (  # noqa: E402
    BENCHMARKS, DTYPE, ROOT, SEED, BenchmarkCase, benchmark_case,
    dense_linear_response, ensure_tree, environment_record,
    formulation_energy, geometry_diagnostics, invariant_metrics,
    mesh_quality, method_unavailable_rows, plot_style, relative_error,
    rss_megabytes, save_figure, state_counts, timed_median, write_csv,
    write_json,
)
from checkerboard_shell import (  # noqa: E402
    BoundaryConditions, CheckerboardMesh, CheckerboardShell,
    IsotropicMaterial, NodalLoad, NewtonSettings,
)
from checkerboard_shell.geometry.core import edge_midpoints  # noqa: E402
from checkerboard_shell.mechanics.energy import make_targets, paper_components  # noqa: E402


torch.set_default_dtype(DTYPE)
torch.manual_seed(SEED)
np.random.seed(SEED)


def log_line(handle, message: str) -> None:
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    text = f"[{stamp}] {message}"
    print(text, flush=True)
    handle.write(text + "\n")
    handle.flush()


def _result_row(case: BenchmarkCase, result: dict, context: dict) -> dict:
    total_states, active_states = state_counts(case)
    row = {
        **context, "method": "CB-Full", "status": "CONVERGED_LINEAR",
        "reference_displacement": case.reference_value,
        "observable": result["observable"],
        "absolute_error": abs(result["observable"]-case.reference_value),
        "relative_error": relative_error(result["observable"], case.reference_value),
        "scalar_state_count": total_states, "active_state_count": active_states,
        "linear_iterations": 1, "newton_iterations": 0,
        "condition_estimate": result["condition_estimate"],
        "tangent_assembly_time_s": result["assembly_time_s"],
        "linear_solve_time_s": result["solve_time_s"],
        "residual_l2": result["residual_l2"], "nullity": result["nullity"],
    }
    row.update(mesh_quality(case.mesh))
    row.update(geometry_diagnostics(result, case))
    positives = result["eigenvalues"][result["eigenvalues"] >
                                      1.0e-10*result["eigenvalues"].abs().max()]
    for index, value in enumerate(positives[:6]):
        row[f"positive_eigenvalue_{index+1}"] = float(value)
    return row


def run_distortion(log) -> list[dict]:
    base = ROOT / "01_mesh_distortion"
    config = {
        "benchmarks": BENCHMARKS,
        "mesh": {"nominal_n": 8, "convergence_n": [4, 6, 8, 10, 12]},
        "smooth_and_random_delta": [0, 0.05, 0.10, 0.20, 0.30, 0.40],
        "nonuniform_spacing_ratio": [1, 2, 4, 8],
        "directional_nx_ny": [[8, 8], [12, 6], [16, 4], [24, 3]],
        "seed": SEED,
        "solution": "geometrically linear reference tangent",
        "baseline_policy": "No unvalidated replacement: absent baselines recorded unavailable.",
    }
    write_json(base / "config" / "experiment.json", config)
    rows: list[dict] = []

    def execute(context: dict, **case_options):
        name = context["benchmark"]
        try:
            case = benchmark_case(name, case_options.pop("nx"), case_options.pop("ny", None),
                                  **case_options)
            quality = mesh_quality(case.mesh)
            if quality["inverted_corner_count"]:
                rows.append({**context, "method": "CB-Full", "status": "INVALID_MESH",
                             **quality, "reason": "non-positive corner scaled Jacobian"})
                return
            result = dense_linear_response(case)
            rows.append(_result_row(case, result, context))
            if (context.get("study") == "distortion"
                    and context.get("distortion_kind") == "smooth"
                    and context.get("delta") == 0.4):
                np.savez_compressed(
                    base / "raw_data" / f"representative_{name}_smooth_040.npz",
                    reference=case.mesh.reference_vertices.numpy(),
                    deformed=result["positions"].numpy(),
                    quads=case.mesh.quads)
        except Exception as error:
            rows.append({**context, "method": "CB-Full", "status": "FAILED",
                         "reason": repr(error)})
            log_line(log, f"Experiment 1 case failed {context}: {error}")
            traceback.print_exc(file=log)

    for benchmark in ("plate", "roof", "m3"):
        for kind in ("smooth", "random"):
            for delta in config["smooth_and_random_delta"]:
                execute({"study": "distortion", "benchmark": benchmark,
                         "distortion_kind": kind, "delta": delta},
                        nx=8, ny=8, distortion=delta,
                        distortion_kind=("none" if delta == 0 else kind))
        for ratio in config["nonuniform_spacing_ratio"]:
            execute({"study": "nonuniform", "benchmark": benchmark,
                     "spacing_ratio": ratio}, nx=8, ny=8,
                    nonuniform_ratio=ratio)
        for nx, ny in config["directional_nx_ny"]:
            execute({"study": "directional", "benchmark": benchmark,
                     "nx": nx, "ny": ny, "directional_ratio": nx/ny},
                    nx=nx, ny=ny)
        for n in config["mesh"]["convergence_n"]:
            execute({"study": "convergence", "benchmark": benchmark, "n": n},
                    nx=n, ny=n)

    contexts = [{key: row.get(key) for key in
                 ("study", "benchmark", "distortion_kind", "delta", "spacing_ratio",
                  "nx", "ny", "directional_ratio", "n") if key in row}
                for row in rows if row.get("method") == "CB-Full"]
    reason_mid = "MidedgeTan is not present in this repository or its declared dependencies."
    reason_mitc = "No verified MITC4 backend is present; an ad-hoc substitute was prohibited."
    rows += method_unavailable_rows("MidedgeTan", contexts, reason_mid)
    rows += method_unavailable_rows("MITC4", contexts, reason_mitc)
    write_csv(base / "raw_data" / "all_cases.csv", rows)
    write_json(base / "raw_data" / "all_cases.json", rows)
    _plot_distortion(rows)
    log_line(log, f"Experiment 1 completed: {len(rows)} rows")
    return rows


def _plot_distortion(rows: list[dict]) -> None:
    plot_style()
    import matplotlib.pyplot as plt

    valid = [r for r in rows if r.get("status") == "CONVERGED_LINEAR"]
    markers = {"plate": "o", "roof": "s", "m3": "^"}
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.7))
    for benchmark in markers:
        subset = [r for r in valid if r.get("study") == "distortion"
                  and r.get("distortion_kind") == "smooth"
                  and r["benchmark"] == benchmark]
        subset.sort(key=lambda r: r["delta"])
        axes[0].plot([r["delta"] for r in subset],
                     [r["relative_error"] for r in subset],
                     marker=markers[benchmark], label=benchmark)
        axes[1].plot([r["minimum_scaled_jacobian"] for r in subset],
                     [r["relative_error"] for r in subset],
                     marker=markers[benchmark], label=benchmark)
    axes[0].set(xlabel=r"distortion $\delta$", ylabel="relative error", yscale="log")
    axes[1].set(xlabel="minimum scaled Jacobian", ylabel="relative error", yscale="log")
    axes[0].legend(frameon=False)
    for ax in axes:
        ax.grid(True, which="both", alpha=0.25)
    save_figure(fig, ROOT / "01_mesh_distortion" / "figures" / "Figure_D1")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(3.5, 2.7))
    for benchmark in markers:
        subset = [r for r in valid if r.get("study") == "nonuniform"
                  and r["benchmark"] == benchmark]
        subset.sort(key=lambda r: r["spacing_ratio"])
        ax.plot([r["spacing_ratio"] for r in subset],
                [r["relative_error"] for r in subset],
                marker=markers[benchmark], label=benchmark)
    ax.set(xlabel="max/min local spacing", ylabel="relative error", xscale="log", yscale="log")
    ax.set_xticks([1, 2, 4, 8], ["1", "2", "4", "8"])
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    save_figure(fig, ROOT / "01_mesh_distortion" / "figures" / "Figure_D2")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(3.5, 2.7))
    for benchmark in markers:
        subset = [r for r in valid if r.get("study") == "convergence"
                  and r["benchmark"] == benchmark]
        subset.sort(key=lambda r: r["scalar_state_count"])
        ax.plot([r["scalar_state_count"] for r in subset],
                [r["relative_error"] for r in subset],
                marker=markers[benchmark], label=f"CB-Full: {benchmark}")
    ax.set(xlabel="total scalar state count", ylabel="relative error", xscale="log", yscale="log")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    save_figure(fig, ROOT / "01_mesh_distortion" / "figures" / "Figure_D3")
    plt.close(fig)

    archives = list((ROOT / "01_mesh_distortion" / "raw_data").glob("representative_*.npz"))
    if archives:
        fig = plt.figure(figsize=(7.1, 2.5))
        for index, archive in enumerate(archives[:3], 1):
            data = np.load(archive)
            ax = fig.add_subplot(1, 3, index, projection="3d")
            ref, deformed = data["reference"], data["deformed"]
            maximum_u = max(float(np.linalg.norm(deformed-ref, axis=-1).max()), 1e-30)
            span = float(np.linalg.norm(np.ptp(ref.reshape(-1, 3), axis=0)))
            visibility_scale = max(1.0, 0.08*span/maximum_u)
            plotted = ref+visibility_scale*(deformed-ref)
            for field, style, alpha in ((ref, "k-", 0.35), (plotted, "C1-", 0.9)):
                for row in field:
                    ax.plot(row[:, 0], row[:, 1], row[:, 2], style, lw=0.45, alpha=alpha)
                for col in field.transpose(1, 0, 2):
                    ax.plot(col[:, 0], col[:, 1], col[:, 2], style, lw=0.45, alpha=alpha)
            ax.set_title(archive.stem.replace("representative_", "")
                         + f"\norange: {visibility_scale:.1f}× displacement")
            ax.set_axis_off()
        save_figure(fig, ROOT / "01_mesh_distortion" / "figures" / "Figure_D4")
        plt.close(fig)


def _large_deformation_state(log):
    mesh = CheckerboardMesh.rectangle(width=2.0, height=0.6, nx=6, ny=4, dtype=DTYPE)
    material = IsotropicMaterial(young=2.0e4, poisson=0.3, thickness=0.08)
    boundary = BoundaryConditions.clamped_boundary(mesh, "left")
    load = NodalLoad.total_force(mesh, (0.0, 0.0, -0.08), mesh.boundaries["right"])
    model = CheckerboardShell(
        mesh, material, boundary, load,
        NewtonSettings(residual_tolerance=1.0e-8, maximum_newton_iterations=30))
    result = model.solve()
    log_line(log, f"Gauge large-deformation cantilever status={result.status}, "
                  f"iterations={result.iterations}, max_u={result.max_displacement:.6e}")
    return mesh, material, result.deformed_vertices, result.status


def run_gauge(log) -> list[dict]:
    base = ROOT / "02_gauge_invariance"
    direction = torch.tensor([1.0, 0.7, -0.4], dtype=DTYPE)
    direction /= torch.linalg.norm(direction)
    multipliers = [0.0, 1e-8, 1e-4, 1e-2, 0.1, 0.5, 1.0]
    write_json(base / "config" / "experiment.json", {
        "alpha_over_characteristic_length": multipliers,
        "direction": direction.tolist(), "states": [
            "undeformed curved M3", "linear M3 deformed", "nonlinear cantilever final"],
    })
    rows = []
    states = []
    m3 = benchmark_case("m3", 6)
    states.append(("undeformed_curved", m3.mesh, m3.material,
                   m3.mesh.reference_vertices, "REFERENCE"))
    m3_result = dense_linear_response(m3)
    states.append(("linear_m3", m3.mesh, m3.material, m3_result["positions"],
                   "CONVERGED_LINEAR"))
    cant_mesh, cant_material, cant_positions, cant_status = _large_deformation_state(log)
    states.append(("nonlinear_cantilever", cant_mesh, cant_material,
                   cant_positions, cant_status))
    raw_plot = []
    for state_name, mesh, material, current, status in states:
        np.savez_compressed(
            base / "raw_data" / f"state_{state_name}.npz",
            reference=mesh.reference_vertices.detach().numpy(),
            current=current.detach().numpy(), quads=mesh.quads)
        extent = torch.linalg.vector_norm(
            mesh.reference_vertices.reshape(-1, 3).max(0).values
            - mesh.reference_vertices.reshape(-1, 3).min(0).values)
        alphas = [value*float(extent) for value in multipliers]
        try:
            state_rows = invariant_metrics(
                mesh.reference_vertices, current, mesh.hx, mesh.hy,
                material, alphas, direction)
            for row, multiplier in zip(state_rows, multipliers):
                rows.append({"state": state_name, "state_status": status,
                             "alpha_over_L": multiplier, **row})
            raw_plot.append((state_name, current.detach(),
                             (current+alphas[-1]*
                              (((torch.arange(mesh.ny+1)[:, None]
                                 + torch.arange(mesh.nx+1)[None, :]) % 2)*-2+1)
                              .to(DTYPE)[..., None]*direction).detach()))
        except Exception as error:
            rows.append({"state": state_name, "state_status": status,
                         "status": "FAILED", "reason": repr(error)})
            log_line(log, f"Gauge state {state_name} failed: {error}")
            traceback.print_exc(file=log)
    write_csv(base / "raw_data" / "invariance_errors.csv", rows)
    write_json(base / "raw_data" / "invariance_errors.json", rows)
    _plot_gauge(rows, raw_plot)
    log_line(log, f"Experiment 2 completed: {len(rows)} rows")
    return rows


def _plot_gauge(rows, raw_plot):
    plot_style()
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(7.1, 2.5))
    if raw_plot:
        name, original, shifted = raw_plot[0]
        for panel, field, title in ((1, original, r"raw $\alpha=0$"),
                                    (2, shifted, r"raw $\alpha=L$")):
            ax = fig.add_subplot(1, 3, panel, projection="3d")
            arr = field.numpy()
            for line in arr:
                ax.plot(*line.T, "C0-", lw=0.5)
            for line in arr.transpose(1, 0, 2):
                ax.plot(*line.T, "C0-", lw=0.5)
            ax.set_title(title); ax.set_axis_off()
        ax = fig.add_subplot(1, 3, 3, projection="3d")
        for field, color, label in ((original, "k", r"$\alpha=0$"),
                                    (shifted, "C3", r"$\alpha=L$")):
            mids = edge_midpoints(field)
            xyz = torch.cat((mids.horizontal.reshape(-1, 3),
                             mids.vertical.reshape(-1, 3))).numpy()
            ax.scatter(*xyz.T, s=2, c=color, label=label, alpha=0.7)
        ax.set_title("physical edge midpoints"); ax.set_axis_off(); ax.legend(frameon=False)
    save_figure(fig, ROOT / "02_gauge_invariance" / "figures" / "Figure_G1")
    plt.close(fig)

    metrics = {
        "physical_midpoint_change": "physical midpoints (abs.)",
        "relative_a_B_change": r"$a_B$ (relative)",
        "relative_b_W_change": r"$b_W$ (relative)",
        "absolute_XM_change": r"$X_M$ (abs.)",
        "absolute_XWrel_change": r"$X_W^{rel}$ (abs.)",
        "relative_energy_change": "internal energy (relative)",
        "relative_residual_change": "internal residual (relative)",
        "relative_tangent_change": "internal tangent (relative)",
    }
    fig, ax = plt.subplots(figsize=(4.8, 3.2))
    valid = [r for r in rows if "alpha_over_L" in r and r["alpha_over_L"] > 0]
    for metric, label in metrics.items():
        subset = [r for r in valid if r["state"] == "linear_m3"]
        ax.plot([r["alpha_over_L"] for r in subset],
                [max(r[metric], 1e-19) for r in subset], marker="o", label=label)
    ax.set(xlabel=r"$\alpha/L$", ylabel="invariance error", xscale="log", yscale="log")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, ncol=2, fontsize=8, loc="upper center",
              bbox_to_anchor=(0.5, -0.23))
    save_figure(fig, ROOT / "02_gauge_invariance" / "figures" / "Figure_G2")
    plt.close(fig)


def _operator_storage_bytes(n: int) -> int:
    return 4*(n+1)*(n+1)*n*n*8


def run_scaling(log) -> list[dict]:
    base = ROOT / "03_scaling"
    sizes = [16, 32, 64, 128]
    variants = ["CB-0", "CB-M", "CB-Full"]
    write_json(base / "config" / "experiment.json", {
        "benchmarks": ["plate", "m3"], "nominal_n": sizes,
        "variants": variants, "timing_repeats": 3, "warmups": 1,
        "explicit_tangent_n": [16],
        "resource_guard": "n>=64 reports analytically audited dense reconstruction storage; n=128 is not allocated",
    })
    rows = []
    for benchmark in ("plate", "m3"):
        for n in sizes:
            predicted_ops = _operator_storage_bytes(n)
            if n >= 64:
                for variant in variants:
                    case = benchmark_case(benchmark, n)
                    total, active = state_counts(case)
                    rows.append({
                        "benchmark": benchmark, "n": n, "method": variant,
                        "status": "RESOURCE_GUARD",
                        "scalar_state_count": total, "active_state_count": active,
                        "dense_boundary_operator_bytes": predicted_ops,
                        "dense_tangent_bytes": active*active*8,
                        "reason": ("Current boundary closure caches four dense reconstruction matrices; "
                                   "allocation/timing was stopped before risking multi-GB memory."),
                    })
                continue
            case = benchmark_case(benchmark, n)
            total, active = state_counts(case)
            try:
                reference = case.mesh.reference_vertices
                targets = make_targets(reference, case.mesh.hx, case.mesh.hy)
                q0 = case.boundary.initial_coordinates(reference)
                vector = torch.randn_like(q0); vector /= torch.linalg.norm(vector)
                for variant in variants:
                    def energy(q):
                        return formulation_energy(
                            case.boundary.positions(reference, q), reference,
                            case.mesh, targets, case.material, variant)[0]

                    def residual_eval():
                        q = q0.detach().clone().requires_grad_(True)
                        value = energy(q)
                        return torch.autograd.grad(value, q)[0]

                    def tangent_action():
                        return torch.autograd.functional.hvp(energy, q0, vector)[1]

                    residual_time, _ = timed_median(residual_eval)
                    hvp_time, _ = timed_median(tangent_action)
                    row = {
                        "benchmark": benchmark, "n": n, "method": variant,
                        "status": "MATRIX_FREE_MEASURED", "scalar_state_count": total,
                        "active_state_count": active, "residual_eval_time_s": residual_time,
                        "tangent_action_time_s": hvp_time,
                        "peak_process_rss_mb": rss_megabytes(),
                        "dense_boundary_operator_bytes": predicted_ops,
                        "dense_tangent_bytes": active*active*8,
                    }
                    if n == 16:
                        tangent_time, tangent = timed_median(
                            lambda: torch.autograd.functional.hessian(
                                energy, q0, vectorize=True), repeats=3, warmups=1)
                        threshold = max(float(tangent.abs().max())*1e-12, 1e-30)
                        nnz = int((tangent.abs() > threshold).sum())
                        density = nnz/tangent.numel()
                        force = case.load[case.boundary.free_mask]
                        # Gauge-safe spectral pseudoinverse is timed once here.
                        started = time.perf_counter()
                        values, vectors_e = torch.linalg.eigh(0.5*(tangent+tangent.T))
                        keep = values > max(float(values.abs().max())*1e-10, 1e-30)
                        solution = vectors_e[:, keep] @ ((vectors_e[:, keep].T@force)/values[keep])
                        solve_time = time.perf_counter()-started
                        displacement = case.boundary.positions(reference, solution)-reference
                        row.update({"explicit_tangent_time_s": tangent_time,
                                    "tangent_nnz": nnz, "tangent_density": density,
                                    "factorization_solve_time_s": solve_time,
                                    "observable": float((case.observable_weights*displacement).sum()),
                                    "relative_error": relative_error(
                                        float((case.observable_weights*displacement).sum()),
                                        case.reference_value)})
                    rows.append(row)
            except Exception as error:
                for variant in variants:
                    rows.append({"benchmark": benchmark, "n": n, "method": variant,
                                 "status": "FAILED", "reason": repr(error)})
                log_line(log, f"Scaling {benchmark} n={n} failed: {error}")
                traceback.print_exc(file=log)

    # Exact rank-one audit at reference state.
    audit = []
    for benchmark in ("plate", "m3"):
        case = benchmark_case(benchmark, 8)
        ref = case.mesh.reference_vertices
        targets = make_targets(ref, case.mesh.hx, case.mesh.hy)
        q0 = case.boundary.initial_coordinates(ref)

        def e0(q):
            return formulation_energy(case.boundary.positions(ref, q), ref, case.mesh,
                                      targets, case.material, "CB-Full", 0.0)[0]
        def e1(q):
            return formulation_energy(case.boundary.positions(ref, q), ref, case.mesh,
                                      targets, case.material, "CB-Full", 1.0)[0]
        k0 = torch.autograd.functional.hessian(e0, q0, vectorize=True)
        k1 = torch.autograd.functional.hessian(e1, q0, vectorize=True)

        def xm(q):
            return paper_components(case.boundary.positions(ref, q), ref,
                                    case.mesh.hx, case.mesh.hy, targets,
                                    case.material)["X_M"]
        variable = q0.detach().clone().requires_grad_(True)
        gx, = torch.autograd.grad(xm(variable), variable)
        coefficient = (case.material.thickness/4.0)*8.0*case.material.lame_mu
        rank_one = 2.0*coefficient*torch.outer(gx, gx)
        difference = k1-k0
        threshold = max(float(difference.abs().max())*1e-12, 1e-30)
        audit.append({
            "benchmark": benchmark, "n": 8,
            "relative_rank_one_error": float(torch.linalg.norm(difference-rank_one)
                                               / torch.linalg.norm(difference)),
            "xm_difference_numerical_rank": int(torch.linalg.matrix_rank(
                difference, atol=threshold, rtol=0.0)),
            "xm_difference_nnz": int((difference.abs() > threshold).sum()),
            "xm_difference_density": float((difference.abs() > threshold).sum()/difference.numel()),
            "explicit_outer_product_bytes": difference.numel()*8,
            "low_rank_vector_bytes": gx.numel()*8,
        })
    contexts = [{"benchmark": benchmark, "n": n}
                for benchmark in ("plate", "m3") for n in sizes]
    rows += method_unavailable_rows(
        "MidedgeTan", contexts,
        "No verified MidedgeTan source or dependency backend exists in this checkout.")
    rows += method_unavailable_rows(
        "MITC4", contexts,
        "No verified MITC4 source or dependency backend exists in this checkout.")
    write_csv(base / "raw_data" / "scaling.csv", rows)
    write_json(base / "raw_data" / "rank_one_audit.json", audit)
    low_rank_rows = _low_rank_benchmark(log)
    write_csv(base / "raw_data" / "low_rank_equivalence.csv", low_rank_rows)
    _plot_scaling(rows)
    log_line(log, f"Experiment 3 completed: {len(rows)} scaling rows")
    return rows


def _low_rank_benchmark(log) -> list[dict]:
    """Exactly equivalent reference-tangent rank-one sparse solve.

    At the stress-free reference X_M=0, D²(c X_M²) = 2c ggᵀ.  The
    Sherman--Morrison solver below stores only g, and the sparse local K.
    One raw gauge representative per unconstrained parity channel is fixed
    identically for both naive and low-rank solves; this is a coordinate
    choice, not an energy modification.
    """
    from scipy.linalg import qr
    from scipy.sparse import bmat, csc_matrix
    from scipy.sparse.linalg import splu

    rows = []
    for benchmark in ("plate", "m3"):
        for n in (8, 16):
            try:
                case = benchmark_case(benchmark, n)
                reference = case.mesh.reference_vertices
                targets = make_targets(reference, case.mesh.hx, case.mesh.hy)
                q0 = case.boundary.initial_coordinates(reference)

                def energy(q, multiplier):
                    return formulation_energy(
                        case.boundary.positions(reference, q), reference,
                        case.mesh, targets, case.material,
                        "CB-Full", multiplier)[0]

                def xm(q):
                    return paper_components(
                        case.boundary.positions(reference, q), reference,
                        case.mesh.hx, case.mesh.hy, targets,
                        case.material)["X_M"]

                k_local = torch.autograd.functional.hessian(
                    lambda q: energy(q, 0.0), q0, vectorize=True).numpy()
                variable = q0.detach().clone().requires_grad_(True)
                gradient, = torch.autograd.grad(xm(variable), variable)
                g = gradient.detach().numpy()
                c = 2.0*(case.material.thickness/4.0)*8.0*case.material.lame_mu
                force = case.load[case.boundary.free_mask].numpy()
                k_all = k_local+c*np.outer(g, g)
                all_values, all_vectors = np.linalg.eigh(k_all)
                null = np.flatnonzero(all_values <= max(all_values[-1]*1e-10, 1e-20))
                gauge = np.zeros(q0.numel(), dtype=bool)
                if len(null):
                    # QR selects independent raw-coordinate representatives
                    # for the actual constrained nullspace.  Both formulations
                    # use the same quotient-space coordinates.
                    _, _, pivots = qr(all_vectors[:, null].T, pivoting=True)
                    gauge[pivots[:len(null)]] = True
                keep = np.flatnonzero(~gauge)
                k0 = k_local[np.ix_(keep, keep)]
                gg = g[keep]
                ff = force[keep]
                k_naive = k0+c*np.outer(gg, gg)
                eig = np.linalg.eigvalsh(k_naive)
                if eig[0] <= max(eig[-1]*1e-11, 1e-20):
                    rows.append({"benchmark": benchmark, "n": n,
                                 "status": "SINGULAR_GAUGE_REDUCTION",
                                 "minimum_eigenvalue": float(eig[0])})
                    continue
                t_naive, q_naive = timed_median(
                    lambda: np.linalg.solve(k_naive, ff), repeats=3, warmups=1)
                threshold = max(float(np.max(np.abs(k0)))*1e-15, 1e-30)
                sparse_local = csc_matrix(
                    np.where(np.abs(k0)>threshold, k0, 0.0))
                # Augmented saddle block is exactly equivalent to the rank-one
                # K_local+c*g*g^T, including when K_local alone is singular.
                augmented = bmat([
                    [sparse_local, csc_matrix(gg[:, None])],
                    [csc_matrix(gg[None, :]), csc_matrix([[-1.0/c]])],
                ], format="csc")
                t_factor, lu = timed_median(
                    lambda: splu(augmented), repeats=3, warmups=1)

                def woodbury():
                    return lu.solve(np.r_[ff, 0.0])[:-1]

                t_update, q_low_rank = timed_median(
                    woodbury, repeats=3, warmups=1)
                full_q = np.zeros(q0.numel()); full_q[keep] = q_low_rank
                naive_q = np.zeros(q0.numel()); naive_q[keep] = q_naive
                displacement = case.boundary.positions(
                    reference, torch.as_tensor(full_q, dtype=DTYPE))-reference
                naive_displacement = case.boundary.positions(
                    reference, torch.as_tensor(naive_q, dtype=DTYPE))-reference
                energy_gap = abs(float(energy(torch.as_tensor(full_q, dtype=DTYPE), 1.0)
                                       -energy(torch.as_tensor(naive_q, dtype=DTYPE), 1.0)))
                rows.append({
                    "benchmark": benchmark, "n": n, "status": "VERIFIED",
                    "states_after_gauge_fix": len(keep), "gauge_coordinates_fixed": int(gauge.sum()),
                    "local_sparse_nnz": sparse_local.nnz,
                    "local_sparse_density": sparse_local.nnz/(len(keep)**2),
                    "naive_full_nnz": int(np.count_nonzero(np.abs(k_naive)>threshold)),
                    "naive_dense_tangent_bytes": k_naive.nbytes,
                    "low_rank_vector_bytes": gg.nbytes,
                    "local_sparse_csc_bytes": (sparse_local.data.nbytes
                                               +sparse_local.indices.nbytes
                                               +sparse_local.indptr.nbytes),
                    "augmented_sparse_nnz": augmented.nnz,
                    "augmented_sparse_csc_bytes": (augmented.data.nbytes
                                                   +augmented.indices.nbytes
                                                   +augmented.indptr.nbytes),
                    "naive_solve_median_s": t_naive,
                    "sparse_factorization_median_s": t_factor,
                    "low_rank_update_solve_median_s": t_update,
                    "relative_displacement_difference": float(
                        np.linalg.norm(q_naive-q_low_rank)
                        /max(np.linalg.norm(q_naive), 1e-300)),
                    "relative_observable_difference": float(
                        abs(((case.observable_weights*(displacement-naive_displacement)).sum()))
                        /max(abs(float((case.observable_weights*naive_displacement).sum())), 1e-300)),
                    "absolute_energy_difference": energy_gap,
                    "relative_linear_residual": float(
                        np.linalg.norm(k_naive@q_low_rank-ff)
                        /max(np.linalg.norm(ff), 1e-300)),
                })
            except Exception as error:
                rows.append({"benchmark": benchmark, "n": n, "status": "FAILED",
                             "reason": repr(error)})
                log_line(log, f"rank-one benchmark {benchmark} n={n} failed: {error}")
                traceback.print_exc(file=log)
    return rows


def _plot_scaling(rows):
    plot_style()
    import matplotlib.pyplot as plt
    styles = {"CB-0": "o-", "CB-M": "s--", "CB-Full": "^-"}
    measured = [r for r in rows if r.get("status") == "MATRIX_FREE_MEASURED"]
    fig, ax = plt.subplots(figsize=(3.5, 2.7))
    for method, style in styles.items():
        subset = [r for r in measured if r["method"] == method and r["benchmark"] == "plate"]
        subset.sort(key=lambda r: r["scalar_state_count"])
        ax.plot([r["scalar_state_count"] for r in subset],
                [r["tangent_action_time_s"] for r in subset], style, label=method)
    ax.set(xlabel="scalar states", ylabel="tangent action time [s]", xscale="log", yscale="log")
    ax.grid(True, which="both", alpha=0.25); ax.legend(frameon=False)
    save_figure(fig, ROOT / "03_scaling" / "figures" / "Figure_C1")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(3.5, 2.7))
    for method, style in styles.items():
        subset = [r for r in rows if r["method"] == method and r["benchmark"] == "plate"]
        subset.sort(key=lambda r: r["scalar_state_count"])
        ax.plot([r["scalar_state_count"] for r in subset],
                [r["dense_boundary_operator_bytes"]/2**20 for r in subset], style, label=method)
    ax.set(xlabel="scalar states", ylabel="boundary-operator storage [MiB]", xscale="log", yscale="log")
    ax.grid(True, which="both", alpha=0.25); ax.legend(frameon=False)
    save_figure(fig, ROOT / "03_scaling" / "figures" / "Figure_C2")
    plt.close(fig)

    explicit = [r for r in measured if isinstance(r.get("tangent_density"), (int, float))
                and r["benchmark"] == "plate"]
    fig, ax = plt.subplots(figsize=(3.5, 2.7))
    ax.bar([r["method"] for r in explicit], [r["tangent_density"] for r in explicit],
           color="0.65", edgecolor="black", hatch=["//", "..", "xx"][:len(explicit)])
    ax.set(ylabel="numerical tangent density", ylim=(0, 1.05))
    ax.grid(True, axis="y", alpha=0.25)
    save_figure(fig, ROOT / "03_scaling" / "figures" / "Figure_C3")
    plt.close(fig)


def run_beta(log) -> list[dict]:
    base = ROOT / "04_beta_sensitivity"
    multipliers = [0, 0.25, 0.5, 1, 2, 4, 8]
    meshes = {"m3": [8, 12], "roof": [8, 12]}
    write_json(base / "config" / "experiment.json", {
        "c_beta": multipliers, "meshes": meshes,
        "definition": "c_beta = beta_M/(8 mu)",
        "solution": "reference-tangent linear response; K(c)=K(0)+c K_XM",
    })
    rows = []
    for benchmark, sizes in meshes.items():
        for n in sizes:
            try:
                case = benchmark_case(benchmark, n)
                ref = case.mesh.reference_vertices
                targets = make_targets(ref, case.mesh.hx, case.mesh.hy)
                q0 = case.boundary.initial_coordinates(ref)

                def energy(q, c):
                    return formulation_energy(
                        case.boundary.positions(ref, q), ref, case.mesh, targets,
                        case.material, "CB-Full", c)[0]
                kbase = torch.autograd.functional.hessian(
                    lambda q: energy(q, 0.0), q0, vectorize=True)
                kfull = torch.autograd.functional.hessian(
                    lambda q: energy(q, 1.0), q0, vectorize=True)
                kxm = kfull-kbase
                force = case.load[case.boundary.free_mask]
                for multiplier in multipliers:
                    tangent = 0.5*(kbase+multiplier*kxm + (kbase+multiplier*kxm).T)
                    values, vectors = torch.linalg.eigh(tangent)
                    scale = max(float(values.abs().max()), 1e-300)
                    keep = values > 1e-10*scale
                    q = vectors[:, keep] @ ((vectors[:, keep].T@force)/values[keep])
                    positions = case.boundary.positions(ref, q)
                    displacement = positions-ref
                    observable = float((case.observable_weights*displacement).sum())
                    _, parts = formulation_energy(
                        positions, ref, case.mesh, targets, case.material,
                        "CB-Full", multiplier)
                    total = float((parts["membrane_b"]+parts["bending"]
                                   + multiplier*parts["membrane_xm"]
                                   + parts["bending_compatibility"]).detach())
                    positive = values[keep]
                    # Unit direction of the rank-one completion is the target z_M channel.
                    evals_xm, evecs_xm = torch.linalg.eigh(0.5*(kxm+kxm.T))
                    zmode = evecs_xm[:, -1]
                    zstiff = float(zmode @ tangent @ zmode)
                    z_added = float(multiplier*(zmode@kxm@zmode))
                    rows.append({
                        "benchmark": benchmark, "n": n, "c_beta": multiplier,
                        "status": "CONVERGED_LINEAR", "observable": observable,
                        "reference": case.reference_value,
                        "relative_error": relative_error(observable, case.reference_value),
                        "E_XM_over_E_int": (multiplier*float(parts["membrane_xm"].detach())/total
                                             if total else 0.0),
                        "smallest_positive_eigenvalue": float(positive[0]) if len(positive) else math.nan,
                        "second_positive_eigenvalue": float(positive[1]) if len(positive) > 1 else math.nan,
                        "z_M_channel_stiffness": zstiff,
                        "z_M_completion_added_stiffness": z_added,
                        "condition_estimate": (float(positive[-1]/positive[0])
                                               if len(positive) else math.inf),
                        "nullity": int((~keep).sum()),
                    })
            except Exception as error:
                rows.append({"benchmark": benchmark, "n": n, "status": "FAILED",
                             "reason": repr(error)})
                log_line(log, f"Beta {benchmark} n={n} failed: {error}")
                traceback.print_exc(file=log)
    write_csv(base / "raw_data" / "beta_sensitivity.csv", rows)
    write_json(base / "raw_data" / "beta_sensitivity.json", rows)
    _plot_beta(rows)
    log_line(log, f"Experiment 4 completed: {len(rows)} rows")
    return rows


def _plot_beta(rows):
    plot_style()
    import matplotlib.pyplot as plt
    markers = {"m3": "o", "roof": "s"}
    valid = [r for r in rows if r.get("status") == "CONVERGED_LINEAR"]
    for figure, field, ylabel, logy in (
        ("Figure_B1", "relative_error", "relative displacement error", True),
        ("Figure_B2", "z_M_completion_added_stiffness", r"completion-added $z_M$ stiffness", True),
        ("Figure_B3", "E_XM_over_E_int", r"$E_{X_M}/E_{int}$", False),
    ):
        fig, ax = plt.subplots(figsize=(3.5, 2.7))
        for benchmark in markers:
            for n in sorted({r["n"] for r in valid if r["benchmark"] == benchmark}):
                subset = [r for r in valid if r["benchmark"] == benchmark and r["n"] == n]
                subset.sort(key=lambda r: r["c_beta"])
                ax.plot([r["c_beta"] for r in subset],
                        [max(r[field], 1e-30) if logy else r[field] for r in subset],
                        marker=markers[benchmark], linestyle="-" if n == max(x["n"] for x in subset) else "--",
                        label=f"{benchmark}, n={n}")
        ax.set(xlabel=r"$c_\beta=\beta_M/(8\mu)$", ylabel=ylabel)
        if logy: ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.25); ax.legend(frameon=False, fontsize=8)
        save_figure(fig, ROOT / "04_beta_sensitivity" / "figures" / figure)
        plt.close(fig)


def run_consistency(log) -> dict:
    """Selected automatic checks required by section I."""
    case = benchmark_case("m3", 6)
    ref = case.mesh.reference_vertices
    targets = make_targets(ref, case.mesh.hx, case.mesh.hy)
    material = case.material

    def full_energy(vertices):
        return paper_components(vertices, ref, case.mesh.hx, case.mesh.hy,
                                targets, material)["total"]
    variable = ref.detach().clone().requires_grad_(True)
    value = full_energy(variable)
    residual, = torch.autograd.grad(value, variable, create_graph=True)
    tangent = torch.autograd.functional.hessian(
        lambda x: full_energy(x.reshape_as(ref)), ref.reshape(-1), vectorize=True)
    translation = torch.tensor([0.37, -0.21, 0.13], dtype=DTYPE)
    angle = 0.731
    axis = torch.tensor([1.0, 2.0, -0.5], dtype=DTYPE); axis /= torch.linalg.norm(axis)
    cross = torch.tensor([[0., -axis[2], axis[1]], [axis[2], 0., -axis[0]],
                          [-axis[1], axis[0], 0.]], dtype=DTYPE)
    rotation = torch.eye(3, dtype=DTYPE)+math.sin(angle)*cross+(1-math.cos(angle))*(cross@cross)
    rigid = ref@rotation.T+translation
    rigid_energy = float(full_energy(rigid).detach())
    direction = torch.randn(ref.numel(), dtype=DTYPE); direction /= torch.linalg.norm(direction)
    eps = 1e-6
    def flat_residual(x):
        x = x.detach().clone().requires_grad_(True)
        val = full_energy(x.reshape_as(ref))
        return torch.autograd.grad(val, x)[0]
    rp = flat_residual(ref.reshape(-1)+eps*direction)
    rm = flat_residual(ref.reshape(-1)-eps*direction)
    fd = (rp-rm)/(2*eps)
    ad = tangent@direction
    check = {
        "reference_energy": float(value.detach()),
        "reference_residual_norm": float(torch.linalg.norm(residual.detach())),
        "rigid_translation_energy": float(full_energy(ref+translation).detach()),
        "finite_rigid_rotation_energy": rigid_energy,
        "tangent_relative_asymmetry": float(torch.linalg.norm(tangent-tangent.T)
                                             / max(float(torch.linalg.norm(tangent)), 1.0)),
        "directional_derivative_relative_error": float(torch.linalg.norm(fd-ad)
                                                        / max(float(torch.linalg.norm(ad)), 1.0)),
        "reference_minimum_frame_det": mesh_quality(case.mesh)["minimum_b_metric_determinant"],
        "reference_normal_flip": False,
    }
    # Reuse the finite gauge check at alpha=L.
    d = torch.tensor([1.0, 0.7, -0.4], dtype=DTYPE); d /= torch.linalg.norm(d)
    gauge = invariant_metrics(ref, ref, case.mesh.hx, case.mesh.hy, material, [1.0], d)[0]
    check["gauge_max_selected_error"] = max(
        value for key, value in gauge.items()
        if key not in ("alpha", "absolute_internal_energy_change",
                       "absolute_internal_residual_change") and math.isfinite(value))
    check["pass"] = bool(
        abs(check["reference_energy"]) < 1e-12
        and check["reference_residual_norm"] < 1e-10
        and abs(check["rigid_translation_energy"]) < 1e-10
        and abs(check["finite_rigid_rotation_energy"]) < 1e-8
        and check["tangent_relative_asymmetry"] < 1e-10
        and check["directional_derivative_relative_error"] < 1e-6
        and check["gauge_max_selected_error"] < 1e-8)
    write_json(ROOT / "consistency_checks.json", check)
    log_line(log, f"Consistency checks pass={check['pass']}: {check}")
    return check


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--only", choices=["distortion", "gauge", "scaling", "beta", "consistency"])
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    torch.set_num_threads(args.threads)
    ensure_tree()
    write_json(ROOT / "environment.json", environment_record())
    log_path = ROOT / "primary_run.log"
    with log_path.open("a") as log:
        log_line(log, f"Primary run start, only={args.only}, threads={args.threads}")
        if args.only in (None, "consistency"): run_consistency(log)
        if args.only in (None, "distortion"): run_distortion(log)
        if args.only in (None, "gauge"): run_gauge(log)
        if args.only in (None, "scaling"): run_scaling(log)
        if args.only in (None, "beta"): run_beta(log)
        log_line(log, "Primary run finished")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
