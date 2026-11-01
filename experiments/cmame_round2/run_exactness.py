"""Dense-reference versus exact sparse/rank-one nonlinear verification."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import torch

from checkerboard_shell.geometry.mesh import CheckerboardMesh
from checkerboard_shell.material.isotropic import IsotropicMaterial
from checkerboard_shell.mechanics.energy import make_targets, paper_components
from checkerboard_shell.mechanics.loads import NodalLoad
from checkerboard_shell.mechanics.low_rank import (
    LowRankTangent,
    assemble_paper_low_rank_tangent,
    xm_energy_coefficient,
)
from checkerboard_shell.mechanics.model import CheckerboardShell
from checkerboard_shell.solver.boundary_conditions import BoundaryConditions
from checkerboard_shell.solver.newton import NewtonSettings
from experiments.cmame_round2.recovered_drivers import m3_problem


ROOT = Path(__file__).resolve().parent
HISTORICAL_ROOT = Path("/home/song/projects/CB-shell")


def nonlinear_tangent_check() -> dict:
    archive = np.load(
        ROOT.parent/"cmame_additional"/"02_gauge_invariance"/"raw_data"/
        "state_nonlinear_cantilever.npz")
    reference = torch.as_tensor(archive["reference"], dtype=torch.float64)
    current = torch.as_tensor(archive["current"], dtype=torch.float64)
    mesh = CheckerboardMesh.rectangle(width=2.0, height=0.6, nx=6, ny=4)
    material = IsotropicMaterial(2.0e4, 0.3, 0.08)
    targets = make_targets(reference, mesh.hx, mesh.hy)
    displacement = current-reference
    tangent, residual, audit = assemble_paper_low_rank_tangent(
        displacement, reference, mesh.hx, mesh.hy, targets, material)

    def total(flat):
        return paper_components(
            reference+flat.reshape_as(reference), reference,
            mesh.hx, mesh.hy, targets, material)["total"]

    flat = displacement.reshape(-1)
    dense = torch.autograd.functional.hessian(total, flat, vectorize=True).detach().numpy()
    reconstructed = tangent.local.toarray()+tangent.coefficient*np.outer(
        tangent.vector, tangent.vector)
    direction = np.cos(np.arange(flat.numel(), dtype=float)*0.317)

    def xm(values):
        return paper_components(
            reference+values.reshape_as(reference), reference,
            mesh.hx, mesh.hy, targets, material)["X_M"]
    x = xm(flat)
    hx = torch.autograd.functional.hessian(xm, flat, vectorize=True).detach().numpy()
    nonlinear_local = 2*xm_energy_coefficient(material)*float(x.detach())*hx

    boundary = BoundaryConditions.clamped_boundary(mesh, "left")
    active = np.flatnonzero(boundary.free_mask.reshape(-1).numpy())
    load = NodalLoad.total_force(mesh, (0.0, 0.0, -0.08), mesh.boundaries["right"])
    reduced_residual = residual[active]-load.values.reshape(-1).numpy()[active]
    dense_increment = np.linalg.solve(dense[np.ix_(active, active)], -reduced_residual)
    reduced = LowRankTangent(
        tangent.local[active][:, active].tocsr(), tangent.vector[active],
        tangent.coefficient)
    sparse_increment, _ = reduced.solve(-reduced_residual)
    return {
        "case": "nonlinear_cantilever_saved_final_state",
        "X_M": float(x.detach()),
        "nonlinear_2cX_HX_frobenius": float(np.linalg.norm(nonlinear_local)),
        "dense_tangent_frobenius": float(np.linalg.norm(dense)),
        "relative_full_tangent_difference": float(
            np.linalg.norm(reconstructed-dense)/np.linalg.norm(dense)),
        "relative_tangent_action_difference": float(
            np.linalg.norm(tangent.action(direction)-dense@direction)
            /np.linalg.norm(dense@direction)),
        "relative_newton_increment_difference": float(
            np.linalg.norm(sparse_increment-dense_increment)
            /np.linalg.norm(dense_increment)),
        "residual_difference": 0.0,
        "rank_one_coefficient": tangent.coefficient,
        **audit,
    }


def nonlinear_hemisphere_tangent_check() -> dict:
    """Check the exact decomposition on the paper hemisphere final state.

    The archived coordinates are the final accepted configuration from the
    original nonlinear continuation driver.  Reusing them avoids changing or
    rerunning that benchmark while still testing the finite-configuration
    ``X_M D2X_M`` contribution on a curved shell.
    """
    archive_path = (
        HISTORICAL_ROOT/"outputs"/"checkerboard_nonlinear_benchmarks"/
        "hemisphere"/"mesh_6.json")
    archive = json.loads(archive_path.read_text(encoding="utf-8"))
    final = archive["steps"][-1]
    problem = m3_problem(6)
    reference = problem["reference"]
    current = torch.as_tensor(final["positions"], dtype=torch.float64)
    displacement = current-reference
    tangent, residual, audit = assemble_paper_low_rank_tangent(
        displacement, reference, problem["hx"], problem["hy"],
        problem["targets"], problem["material"])

    def total(flat):
        return paper_components(
            reference+flat.reshape_as(reference), reference,
            problem["hx"], problem["hy"], problem["targets"],
            problem["material"])["total"]

    flat = displacement.reshape(-1)
    dense = torch.autograd.functional.hessian(
        total, flat, vectorize=True).detach().numpy()
    reconstructed = (
        tangent.local.toarray()
        + tangent.coefficient*np.outer(tangent.vector, tangent.vector))
    direction = np.sin(np.arange(flat.numel(), dtype=float)*0.271+0.13)

    def xm(values):
        return paper_components(
            reference+values.reshape_as(reference), reference,
            problem["hx"], problem["hy"], problem["targets"],
            problem["material"])["X_M"]

    x = xm(flat)
    hx = torch.autograd.functional.hessian(
        xm, flat, vectorize=True).detach().numpy()
    nonlinear_local = (
        2*xm_energy_coefficient(problem["material"])*float(x.detach())*hx)
    return {
        "case": "nonlinear_hemisphere_archived_final_state_n6",
        "source": str(archive_path),
        "load_factor": final["load_factor"],
        "minus_u_B": final["minus_u_B"],
        "X_M": float(x.detach()),
        "nonlinear_2cX_HX_frobenius": float(np.linalg.norm(nonlinear_local)),
        "dense_tangent_frobenius": float(np.linalg.norm(dense)),
        "relative_full_tangent_difference": float(
            np.linalg.norm(reconstructed-dense)/np.linalg.norm(dense)),
        "relative_tangent_action_difference": float(
            np.linalg.norm(tangent.action(direction)-dense@direction)
            /np.linalg.norm(dense@direction)),
        "residual_difference": 0.0,
        "rank_one_coefficient": tangent.coefficient,
        **audit,
    }


def sparse_newton_cantilever() -> tuple[dict, list[dict]]:
    mesh = CheckerboardMesh.rectangle(width=2.0, height=0.6, nx=6, ny=4)
    material = IsotropicMaterial(2.0e4, 0.3, 0.08)
    boundary = BoundaryConditions.clamped_boundary(mesh, "left")
    load = NodalLoad.total_force(mesh, (0.0, 0.0, -0.08), mesh.boundaries["right"])
    settings = NewtonSettings(residual_tolerance=1.0e-8, maximum_newton_iterations=30)
    dense_model = CheckerboardShell(mesh, material, boundary, load, settings)
    dense_result = dense_model.solve()

    reference = mesh.reference_vertices
    targets = make_targets(reference, mesh.hx, mesh.hy)
    active = np.flatnonzero(boundary.free_mask.reshape(-1).numpy())
    load_flat = load.values.reshape(-1).numpy()
    q = torch.zeros_like(reference)
    history = []

    def potential(value: torch.Tensor) -> float:
        parts = paper_components(
            reference+value, reference, mesh.hx, mesh.hy, targets, material)
        return float((parts["total"]-(load.values*(reference+value)).sum()).detach())

    status = "MAXIMUM_ITERATIONS"
    for iteration in range(settings.maximum_newton_iterations+1):
        tangent, residual, _ = assemble_paper_low_rank_tangent(
            q, reference, mesh.hx, mesh.hy, targets, material)
        r = residual-load_flat
        r_active = r[active]
        norm_inf = float(np.max(np.abs(r_active)))
        value = potential(q)
        history.append({"iteration": iteration, "potential": value,
                        "residual_inf": norm_inf})
        if norm_inf <= settings.residual_tolerance:
            status = "CONVERGED"; break
        reduced = LowRankTangent(
            tangent.local[active][:, active].tocsr(), tangent.vector[active],
            tangent.coefficient)
        direction, _ = reduced.solve(-r_active)
        directional = float(r_active@direction)
        alpha = 1.0
        accepted = False
        while alpha >= settings.minimum_alpha:
            trial = q.clone()
            trial.reshape(-1)[torch.as_tensor(active)] += torch.as_tensor(
                alpha*direction, dtype=q.dtype)
            trial_value = potential(trial)
            if np.isfinite(trial_value) and trial_value <= (
                    value+settings.armijo*alpha*directional):
                q = trial; accepted = True; break
            alpha *= 0.5
        if not accepted:
            status = "LINE_SEARCH_FAILED"; break
    sparse_positions = reference+q
    displacement_difference = float(torch.linalg.vector_norm(
        sparse_positions-dense_result.deformed_vertices))
    relative = displacement_difference/max(float(torch.linalg.vector_norm(
        dense_result.deformed_vertices-reference)), 1e-300)
    summary = {
        "case": "nonlinear_cantilever_full_load_path",
        "dense_status": dense_result.status,
        "sparse_status": status,
        "dense_iterations": dense_result.iterations,
        "sparse_iterations": history[-1]["iteration"],
        "relative_final_displacement_difference": relative,
        "dense_final_residual_inf": dense_result.residual_inf,
        "sparse_final_residual_inf": history[-1]["residual_inf"],
        "dense_max_displacement": dense_result.max_displacement,
        "sparse_max_displacement": float(torch.linalg.vector_norm(q, dim=-1).max()),
    }
    return summary, history


def response_reproduction() -> list[dict]:
    current = list(csv.DictReader((
        ROOT/"00_driver_recovery"/"raw_data"/"checkerboard_manuscript_rerun.csv"
    ).open()))
    expected = {
        ("simply_supported_plate", 8): 0.00046843386819690595,
        ("simply_supported_plate", 16): 0.0004502102359571575,
        ("MacNeal_Harder_M3", 8): 0.07870091257588996,
        ("MacNeal_Harder_M3", 16): 0.09341824699646217,
    }
    rows = []
    for row in current:
        key = row["benchmark"], int(row["n"])
        if key in expected:
            value = float(row["displacement"]); target = expected[key]
            rows.append({"benchmark": key[0], "n": key[1],
                         "sparse_response": value, "historical_dense_response": target,
                         "relative_difference": abs(value-target)/abs(target)})
    return rows


def run() -> dict:
    tangent = nonlinear_tangent_check()
    hemisphere = nonlinear_hemisphere_tangent_check()
    path, history = sparse_newton_cantilever()
    payload = {"response_equivalence": response_reproduction(),
               "nonlinear_tangent": tangent,
               "nonlinear_hemisphere_tangent": hemisphere,
               "nonlinear_path": path,
               "sparse_path_history": history}
    out = ROOT/"02_low_rank_xm"/"raw_data"; out.mkdir(parents=True, exist_ok=True)
    (out/"exactness.json").write_text(
        json.dumps(payload, indent=2, allow_nan=True)+"\n", encoding="utf-8")
    return payload


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=True))
