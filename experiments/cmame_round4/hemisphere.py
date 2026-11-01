"""Large-deformation hemisphere under four alternating radial forces.

Ported from ``CB-shell/experiments/checkerboard_nonlinear_benchmarks/
hemisphere.py`` at tag ``v1.0-arxiv``.  The quarter-symmetry reduction,
load potential, and reference curves are unchanged; the energy is the
current ``checkerboard_shell`` production energy.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import scipy.linalg
import torch

from checkerboard_shell.mechanics.energy import make_targets
from experiments.cmame_round4.common import energy_parts, write_csv, write_json
from experiments.cmame_round4.nonlinear_solver import (
    SolverSettings, adaptive_continuation, interpolate_curve,
)


YOUNG = 6.825e7
POISSON = 0.3
THICKNESS = 0.04
RADIUS = 10.0
MAXIMUM_LOAD = 100.0

REFERENCE_LOADS = np.linspace(0.05, 1.0, 20)
REFERENCE_A_DIGITIZED = np.array([
    .20, .42, .65, .85, 1.03, 1.20, 1.35, 1.50, 1.63, 1.76,
    1.86, 1.97, 2.07, 2.16, 2.25, 2.33, 2.41, 2.48, 2.55, 2.61])
REFERENCE_B_DIGITIZED = np.array([
    .23, .48, .71, .94, 1.17, 1.40, 1.61, 1.82, 2.03, 2.22,
    2.41, 2.59, 2.76, 2.93, 3.10, 3.26, 3.41, 3.56, 3.70, 3.8735])


def reference_grid(n, dtype=torch.float64):
    theta0 = np.deg2rad(18.0)
    theta = torch.linspace(theta0, np.pi/2.0, n + 1, dtype=dtype)
    phi = torch.linspace(0.0, np.pi/2.0, n + 1, dtype=dtype)
    tt, pp = torch.meshgrid(theta, phi, indexing="ij")
    q = RADIUS*torch.stack((
        torch.sin(tt)*torch.cos(pp), torch.sin(tt)*torch.sin(pp),
        torch.cos(tt)), -1)
    return q, float((np.pi/2.0)/n), float((np.pi/2.0-theta0)/n)


def _vertical_midpoint_operator(n: int, i: int, j: int) -> np.ndarray:
    count = 3*(n + 1)**2
    operator = np.zeros((3, count), dtype=float)
    for row in (j, j + 1):
        vertex = row*(n + 1) + i
        for component in range(3):
            operator[component, 3*vertex + component] += 0.5
    return operator


def _equator_point_operator(n: int, i: int) -> np.ndarray:
    if n < 3:
        raise ValueError("hemisphere benchmark requires at least three cells")
    coefficients = (15.0/8.0, -5.0/4.0, 3.0/8.0)
    return sum(
        coefficient*_vertical_midpoint_operator(n, i, n - 1 - offset)
        for offset, coefficient in enumerate(coefficients))


def checkerboard_symmetry_reduction(reference: torch.Tensor):
    """Gauge-native quarter-hemisphere displacement reduction."""
    n = reference.shape[0] - 1
    count = reference.numel()
    constraints = []

    def component_row(operator, component):
        return operator[component].copy()

    for j in range(n):
        left = [_vertical_midpoint_operator(n, i, j) for i in (0, 1, 2)]
        right = [_vertical_midpoint_operator(n, i, j)
                 for i in (n, n - 1, n - 2)]
        constraints.append(component_row(left[0], 1))
        constraints.append(component_row(right[0], 0))
        for component in (0, 2):
            constraints.append(component_row(
                -3.0*left[0] + 4.0*left[1] - left[2], component))
        for component in (1, 2):
            constraints.append(component_row(
                3.0*right[0] - 4.0*right[1] + right[2], component))

    parity = np.fromfunction(
        lambda j, i: (-1.0)**(i + j), (n + 1, n + 1), dtype=int).reshape(-1)
    for component in range(3):
        row = np.zeros(count, dtype=float)
        row[component::3] = parity/parity.size
        constraints.append(row)

    mean_z = np.zeros(count, dtype=float)
    for j in range(n):
        for i in range(n + 1):
            mean_z += component_row(_vertical_midpoint_operator(n, i, j), 2)
    constraints.append(mean_z/(n*(n + 1)))

    constraint = np.stack(constraints)
    transform = scipy.linalg.null_space(constraint, rcond=1.0e-12)
    residual = np.linalg.norm(constraint @ transform, ord=np.inf)
    rank = int(np.linalg.matrix_rank(constraint, tol=1.0e-11))
    return {
        "constraint": constraint,
        "transform": torch.as_tensor(transform, dtype=reference.dtype,
                                     device=reference.device),
        "constraint_rank": rank,
        "constraint_residual_inf": float(residual),
        "point_a": torch.as_tensor(_equator_point_operator(n, 0),
                                   dtype=reference.dtype,
                                   device=reference.device),
        "point_b": torch.as_tensor(_equator_point_operator(n, n),
                                   dtype=reference.dtype,
                                   device=reference.device),
    }


class HemisphereProblem:
    def __init__(self, n: int):
        self.n = n
        self.reference, self.hx, self.hy = reference_grid(n)
        self.targets = make_targets(self.reference, self.hx, self.hy)
        self.reduction = checkerboard_symmetry_reduction(self.reference)
        self.transform = self.reduction["transform"]
        self.base = self.reference.reshape(-1).detach().clone()
        self.point_a = self.reduction["point_a"]
        self.point_b = self.reduction["point_b"]
        self.reference_a = self.point_a @ self.base
        self.reference_b = self.point_b @ self.base
        self.external = (0.5*MAXIMUM_LOAD*self.point_a[0]
                         - 0.5*MAXIMUM_LOAD*self.point_b[1])
        self.initial = torch.zeros(self.transform.shape[1],
                                   dtype=self.reference.dtype)

    def current(self, values: torch.Tensor) -> torch.Tensor:
        return (self.base+self.transform@values).reshape_as(self.reference)

    def parts(self, values: torch.Tensor) -> dict:
        current = self.current(values)
        parts = energy_parts(
            current, self.reference, self.hx, self.hy, self.targets,
            YOUNG, POISSON, THICKNESS)
        parts["_vertices"] = current
        return parts

    def potential(self, values: torch.Tensor, load_factor: float):
        current = self.current(values)
        internal = energy_parts(
            current, self.reference, self.hx, self.hy, self.targets,
            YOUNG, POISSON, THICKNESS)["total"]
        return internal-float(load_factor)*(
            self.external@current.reshape(-1))

    def observables(self, values: torch.Tensor):
        current = self.current(values).reshape(-1)
        a = self.point_a@current
        b = self.point_b@current
        return float(a[0]-self.reference_a[0]), float(
            self.reference_b[1]-b[1])


def _energy_fields(parts: dict) -> dict:
    membrane_b = float(parts["membrane_b"].detach())
    xm = float(parts["membrane_xm"].detach())
    bending = float(parts["bending"].detach())
    xw = float(parts["bending_compatibility"].detach())
    total = float(parts["total"].detach())
    membrane = membrane_b + xm
    return {
        "membrane_B_energy": membrane_b,
        "XM_energy": xm,
        "membrane_energy": membrane,
        "bending_BW_energy": bending,
        "XW_energy": xw,
        "bending_energy": bending+xw,
        "internal_energy": total,
        "XM_fraction_membrane": xm/max(abs(membrane), 1e-30),
        "XM_fraction_internal": xm/max(abs(total), 1e-30),
        "XW_fraction_bending": xw/max(abs(bending+xw), 1e-30),
        "XW_fraction_internal": xw/max(abs(total), 1e-30),
    }


def solve(n: int, output: Path | None = None,
          run_small_gate: bool = False) -> dict:
    problem = HemisphereProblem(n)
    del run_small_gate
    settings = SolverSettings(initial_step=0.01, maximum_step=0.02,
                              minimum_step=0.00125,
                              residual_tolerance=1e-4,
                              maximum_newton_iterations=30,
                              eigenvalue_relative_floor=1e-11)

    def on_accept(load, values, result, step):
        ua, minus_ub = problem.observables(values)
        parts = problem.parts(values)
        row = {
            "load_factor": float(load),
            "load": float(load*MAXIMUM_LOAD),
            "u_A": ua, "minus_u_B": minus_ub,
            "newton_iterations": result["iterations"],
            "minimum_line_search_alpha": result["minimum_alpha"],
            "accepted_step_size": float(step),
            "solver_status": result["status"],
            "residual_inf": result["history"][-1]["residual_inf"],
        }
        row.update(_energy_fields(parts))
        return row

    continuation = adaptive_continuation(
        lambda load, _accepted: lambda x: problem.potential(x, load),
        problem.initial, on_accept, settings,
        mandatory_loads=tuple(REFERENCE_LOADS))
    rows = continuation["accepted"]
    completed = continuation["status"] == "CONVERGED"
    curve = interpolate_curve(rows, REFERENCE_LOADS, ("u_A", "minus_u_B"))
    difference = np.concatenate((
        curve["u_A"]-REFERENCE_A_DIGITIZED,
        curve["minus_u_B"]-REFERENCE_B_DIGITIZED))
    reference = np.concatenate((
        REFERENCE_A_DIGITIZED, REFERENCE_B_DIGITIZED))
    final_a = rows[-1]["u_A"] if rows else None
    final_b = rows[-1]["minus_u_B"] if rows else None
    summary = {
        "case": "hemisphere_alternating_radial_forces",
        "classification": ("VALIDATED_WITH_FINITE_MESH_ERROR"
                           if completed else "SOLVER_LIMITED"),
        "mesh_cells": n,
        "vertices": (n+1)**2,
        "scalar_dof": 3*(n+1)**2,
        "active_dof": problem.initial.numel(),
        "parameters": {"R": RADIUS, "hole_angle_degrees": 18.0,
                       "E": YOUNG, "nu": POISSON, "t": THICKNESS,
                       "Pmax": MAXIMUM_LOAD},
        "formulation_fingerprint": {
            "first_form": "a_B=T_B^T T_B",
            "membrane_energy": "t/4*[sum_B A_B^0 Q(Delta a_B)+8*mu*X_M^2]",
            "second_form": "accepted symmetric b_W",
            "bending_compatibility": "accepted objective/reference-relative X_W_rel",
            "state": "positions only modulo exact checkerboard gauge",
            "energy_source": "checkerboard_shell.mechanics.energy.paper_components",
        },
        "boundary_condition": (
            "physical-midpoint quarter reflection symmetry; gauge slice separate"),
        "load_potential": (
            "dead radial point loads; half resultant at each quarter-model "
            "symmetry point"),
        "continuation_status": continuation["status"],
        "last_converged_load": continuation["last_converged_load"],
        "final_u_A": final_a,
        "final_minus_u_B": final_b,
        "reference_final_minus_u_B": 3.8796,
        "final_A_error_vs_digitized": (
            abs(final_a-REFERENCE_A_DIGITIZED[-1]) /
            REFERENCE_A_DIGITIZED[-1] if final_a else None),
        "final_B_error": (abs(final_b-3.8796)/3.8796 if final_b else None),
        "path_RMS_normalized_error": (
            float(np.linalg.norm(difference)/np.linalg.norm(reference))
            if completed else None),
        "path_max_normalized_error": (
            float(np.max(np.abs(difference))/max(reference))
            if completed else None),
        "reference_curve_note": (
            "published p1 h=0.25 Fig. 3.12 digitization, uncertainty about "
            "0.03; exact converged final B from Table 3.6"),
        "steps": rows,
        "rejected_steps": continuation["rejected"],
        "arc_length_used": False,
    }
    destination = (output or Path(__file__).resolve().parent
                   / "raw_data" / "hemisphere")
    write_json(destination/f"mesh_{n}.json", summary)
    write_csv(destination/f"mesh_{n}_raw_steps.csv", rows)
    return summary


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--meshes", nargs="+", type=int,
                        default=(6, 8, 12, 16))
    arguments = parser.parse_args()
    payload = [solve(n) for n in arguments.meshes]
    print([(p["mesh_cells"], p["classification"],
            p.get("final_B_error"), p.get("path_RMS_normalized_error"))
           for p in payload])
