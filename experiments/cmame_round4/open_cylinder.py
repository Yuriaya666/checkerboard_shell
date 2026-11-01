"""Open-ended cylindrical-shell pull-out with the frozen paper energy.

The adopted problem is the one-eighth Sze--Liu--Lo benchmark reproduced by
Duong--Roohbakhshan--Sauer (CMAME 316, 2017, Sec. 5.2.6).  Only the common
adaptive load-control Newton solver is used; this module contains no arc-
length or formulation-specific stabilization.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import scipy.linalg
import torch

from checkerboard_shell.mechanics.energy import make_targets
from experiments.cmame_round4.common import energy_parts, write_csv, write_json
from experiments.cmame_round4.nonlinear_solver import (
    SolverSettings, adaptive_continuation, solve_load_level)


FORMULATION_FINGERPRINT = {
    "first_form": "a_B=T_B^T T_B",
    "membrane_energy": "t/4*[sum_B A_B^0 Q(Delta a_B)+8*mu*X_M^2]",
    "second_form": "accepted symmetric b_W",
    "bending_compatibility": "accepted objective/reference-relative X_W_rel",
    "state": "positions only modulo exact checkerboard gauge",
    "energy_source": "checkerboard_shell.mechanics.energy.paper_components",
}
OUTPUT = Path(__file__).resolve().parent/"raw_data"/"open_cylinder"


E = 10.5e6
NU = 0.3125
THICKNESS = 0.094
LENGTH = 10.35
RADIUS = 4.953
PMAX = 40_000.0                 # magnitude of either full-cylinder force
REDUCED_FORCE = PMAX/4.0        # one-eighth model share at A

# Vector extraction of the three quartic curves in Duong et al., Fig. 12c.
# These are published high-order reproductions of Sze et al., not an exact
# table.  The independent Sze marker coordinates are stored below.
REFERENCE_LOADS = np.linspace(0.0, 1.0, 41)
REFERENCE_A = np.array([
    0., .818137, 1.257370, 1.523818, 1.702580, 1.831575, 1.929416,
    2.006654, 2.069857, 2.122909, 2.168594, 2.208607, 2.244339,
    2.276488, 2.305885, 2.333054, 2.358333, 2.382121, 2.405011,
    2.427505, 2.451101, 2.519984, 2.551805, 2.573455, 2.590758,
    2.605407, 2.618180, 2.629586, 2.639896, 2.649385, 2.658177,
    2.666273, 2.673971, 2.681162, 2.687971, 2.694484, 2.700613,
    2.706536, 2.712141, 2.717577, 2.722803])
REFERENCE_B = np.array([
    0., .864514, 1.471353, 1.901330, 2.217345, 2.456721, 2.641455,
    2.787172, 2.905219, 3.003956, 3.088560, 3.162712, 3.229001,
    3.289218, 3.344822, 3.396774, 3.446338, 3.494509, 3.542782,
    3.594242, 3.659045, 4.068233, 4.173779, 4.233105, 4.276304,
    4.310694, 4.339544, 4.364535, 4.386591, 4.406547, 4.424630,
    4.441236, 4.456696, 4.471104, 4.484583, 4.497373, 4.509383,
    4.520888, 4.531794, 4.542234, 4.552304])
REFERENCE_C = np.array([
    0., .872576, 1.493946, 1.946615, 2.293881, 2.570282, 2.791444,
    2.966524, 3.105573, 3.217748, 3.309818, 3.386758, 3.451854,
    3.507495, 3.555425, 3.596725, 3.632352, 3.662704, 3.687882,
    3.707484, 3.718231, 3.580300, 3.520506, 3.483975, 3.456486,
    3.434231, 3.415361, 3.398945, 3.384388, 3.371349, 3.359587,
    3.348705, 3.338702, 3.329470, 3.320786, 3.312682, 3.305026,
    3.297768, 3.290910, 3.284314, 3.278093])

# Sze reference markers extracted from the same vector figure.  The marker at
# zero is omitted from metrics and near-axis clipping is retained explicitly.
SZE_MARKERS = {
    "w_A": np.array([
        [.023291565, .819747190], [.076094541, 1.522344982],
        [.150641705, 1.939882553], [.250037927, 2.174081817],
        [.349424032, 2.306459643], [.448820251, 2.418532896],
        [.526472284, 2.540658907], [.604114203, 2.622175776],
        [.697310801, 2.662884448], [.802906640, 2.713844933],
        [.902302862, 2.734149498], [1.0, 2.764705883]]),
    "minus_u_B": np.array([
        [.023295610, .860455860], [.072993720, 1.909346074],
        [.150641706, 2.642530110], [.250037926, 3.080382207],
        [.352528900, 3.334975616], [.448824297, 3.538648355],
        [.517151613, 3.844182346], [.551321339, 4.159828806],
        [.653818380, 4.332984973], [.750107710, 4.414402314],
        [.852604751, 4.485677318], [.951994903, 4.526445710]]),
    "minus_u_C": np.array([
        [.023286507, .860455860], [.072984616, 1.950084605],
        [.150641707, 2.795262272], [.253147852, 3.314621283],
        [.299720865, 3.457201154], [.349418974, 3.559022599],
        [.405336933, 3.630287649], [.455035043, 3.681248137],
        [.650702387, 3.396088385], [.753223702, 3.345177668],
        [.899192938, 3.284064899], [1.0, 3.243356228]]),
}

CASE_OUTPUT = OUTPUT/"open_cylinder_pullout"


def reference_grid(na: int, nc: int, dtype=torch.float64):
    axial = torch.linspace(0.0, LENGTH/2.0, na+1, dtype=dtype)
    theta = torch.linspace(0.0, np.pi/2.0, nc+1, dtype=dtype)
    xx, tt = torch.meshgrid(axial, theta, indexing="ij")
    q = torch.stack((xx, RADIUS*torch.cos(tt), RADIUS*torch.sin(tt)), -1)
    return q, float((LENGTH/2.0)/na), float((np.pi/2.0)/nc)


def _axial_edge(na: int, nc: int, i: int, j: int) -> np.ndarray:
    """Midpoint of the raw edge from axial row j to j+1 at column i."""
    operator = np.zeros((3, 3*(na+1)*(nc+1)))
    for row in (j, j+1):
        vertex = row*(nc+1)+i
        for component in range(3):
            operator[component, 3*vertex+component] += .5
    return operator


def _circumferential_edge(na: int, nc: int, i: int, j: int) -> np.ndarray:
    """Midpoint of the raw edge from circumferential column i to i+1."""
    operator = np.zeros((3, 3*(na+1)*(nc+1)))
    for column in (i, i+1):
        vertex = j*(nc+1)+column
        for component in range(3):
            operator[component, 3*vertex+component] += .5
    return operator


def _axial_endpoint(na: int, nc: int, column: int, at_free_end: bool) -> np.ndarray:
    coefficients = (15/8, -5/4, 3/8)
    rows = ((na-1, na-2, na-3) if at_free_end else (0, 1, 2))
    return sum(coefficient*_axial_edge(na, nc, column, row)
               for coefficient, row in zip(coefficients, rows))


def checkerboard_symmetry_reduction(reference: torch.Tensor) -> dict:
    """Three reflection boundaries plus a distinct exact-gauge slice."""
    na, nc = reference.shape[0]-1, reference.shape[1]-1
    if min(na, nc) < 3:
        raise ValueError("complete B/W stencils require at least 3x3 cells")
    count = reference.numel()
    physical = []

    def component(operator, value):
        return operator[value].copy()

    # theta=0 (z=0) and theta=pi/2 (y=0) reflection planes.
    for j in range(na):
        lower = [_axial_edge(na, nc, i, j) for i in (0, 1, 2)]
        upper = [_axial_edge(na, nc, i, j) for i in (nc, nc-1, nc-2)]
        physical.append(component(lower[0], 2))
        for value in (0, 1):
            physical.append(component(-3*lower[0]+4*lower[1]-lower[2], value))
        physical.append(component(upper[0], 1))
        for value in (0, 2):
            physical.append(component(3*upper[0]-4*upper[1]+upper[2], value))

    # x=0 midspan reflection.  The x=L/2 edge is physically free.
    for i in range(nc):
        edge = [_circumferential_edge(na, nc, i, j) for j in (0, 1, 2)]
        physical.append(component(edge[0], 0))
        for value in (1, 2):
            physical.append(component(-3*edge[0]+4*edge[1]-edge[2], value))

    physical_matrix = np.stack(physical)
    parity = np.fromfunction(
        lambda j, i: (-1.)**(i+j), (na+1, nc+1), dtype=int).reshape(-1)
    gauge = []
    constraints = list(physical)
    for value in range(3):
        row = np.zeros(count)
        row[value::3] = parity/parity.size
        gauge.append(row)
        constraints.append(row)
    gauge = np.column_stack(gauge)
    matrix = np.stack(constraints)
    transform = scipy.linalg.null_space(matrix, rcond=1e-12)
    point_a = _axial_endpoint(na, nc, nc, False)
    point_b = _axial_endpoint(na, nc, 0, False)
    point_c = _axial_endpoint(na, nc, 0, True)
    return {
        "transform": torch.as_tensor(transform, dtype=reference.dtype),
        "point_a": torch.as_tensor(point_a, dtype=reference.dtype),
        "point_b": torch.as_tensor(point_b, dtype=reference.dtype),
        "point_c": torch.as_tensor(point_c, dtype=reference.dtype),
        "constraint_rank": int(np.linalg.matrix_rank(matrix, tol=1e-11)),
        "active_dof": int(transform.shape[1]),
        "physical_constraint_count": len(physical),
        "constraint_residual_inf": float(np.linalg.norm(matrix@transform, ord=np.inf)),
        "gauge_physical_residual_inf": float(
            np.linalg.norm(physical_matrix@gauge, ord=np.inf)),
        "point_gauge_residual_inf": float(max(
            np.linalg.norm(point_a@gauge, ord=np.inf),
            np.linalg.norm(point_b@gauge, ord=np.inf),
            np.linalg.norm(point_c@gauge, ord=np.inf))),
    }


class Problem:
    def __init__(self, na: int, nc: int):
        self.na, self.nc = na, nc
        self.reference, self.hx, self.hy = reference_grid(na, nc)
        self.targets = make_targets(self.reference, self.hx, self.hy)
        self.reduction = checkerboard_symmetry_reduction(self.reference)
        self.transform = self.reduction["transform"]
        self.base = self.reference.reshape(-1).clone()
        self.point_a = self.reduction["point_a"]
        self.point_b = self.reduction["point_b"]
        self.point_c = self.reduction["point_c"]
        self.rest_a = self.point_a@self.base
        self.rest_b = self.point_b@self.base
        self.rest_c = self.point_c@self.base
        self.external = REDUCED_FORCE*self.point_a[2]
        self.initial = torch.zeros(self.transform.shape[1], dtype=torch.float64)

    def current(self, values):
        return (self.base+self.transform@values).reshape_as(self.reference)

    def parts(self, values):
        q = self.current(values)
        parts = energy_parts(
            q, self.reference, self.hx, self.hy, self.targets,
            E, NU, THICKNESS)
        parts["_vertices"] = q
        return parts

    def potential(self, values, load):
        parts = self.parts(values)
        factor = torch.as_tensor(load, dtype=values.dtype, device=values.device)
        return parts["total"]-factor*(self.external@parts["_vertices"].reshape(-1))

    def observables(self, values):
        q = self.current(values).reshape(-1)
        a, b, c = self.point_a@q, self.point_b@q, self.point_c@q
        return (float(a[2]-self.rest_a[2]),
                float(self.rest_b[1]-b[1]),
                float(self.rest_c[1]-c[1]))


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


def _path_metrics(rows: list[dict]) -> dict:
    loads = np.array([0.]+[row["load_factor"] for row in rows])
    values = {
        "w_A": np.array([0.]+[row["w_A"] for row in rows]),
        "minus_u_B": np.array([0.]+[row["minus_u_B"] for row in rows]),
        "minus_u_C": np.array([0.]+[row["minus_u_C"] for row in rows]),
    }
    references = {"w_A": REFERENCE_A, "minus_u_B": REFERENCE_B,
                  "minus_u_C": REFERENCE_C}
    result, all_delta, all_reference = {}, [], []
    for key in values:
        interpolated = np.interp(REFERENCE_LOADS, loads, values[key])
        delta = interpolated-references[key]
        result[key] = {
            "final_relative_error": float(abs(delta[-1])/abs(references[key][-1])),
            "RMS_normalized_path_error": float(
                np.linalg.norm(delta)/np.linalg.norm(references[key])),
            "max_normalized_path_error": float(
                np.max(np.abs(delta))/np.max(np.abs(references[key]))),
        }
        marker = SZE_MARKERS[key]
        marker_delta = np.interp(marker[:, 0], loads, values[key])-marker[:, 1]
        result[key]["Sze_marker_RMS_normalized_error"] = float(
            np.linalg.norm(marker_delta)/np.linalg.norm(marker[:, 1]))
        result[key]["Sze_marker_max_normalized_error"] = float(
            np.max(np.abs(marker_delta))/np.max(np.abs(marker[:, 1])))
        all_delta.append(delta); all_reference.append(references[key])
    delta = np.concatenate(all_delta); reference = np.concatenate(all_reference)
    result["combined"] = {
        "RMS_normalized_path_error": float(np.linalg.norm(delta)/np.linalg.norm(reference)),
        "max_normalized_path_error": float(
            np.max(np.abs(delta))/np.max(np.abs(reference))),
    }
    return result


def solve(na: int, nc: int, *, run_small_gate: bool = False,
          output: Path | None = None) -> dict:
    problem = Problem(na, nc)
    gate = {"passed": None, "rows": []}
    del run_small_gate  # the archived small-load gate is not rerun here
    destination = output or CASE_OUTPUT
    stem = f"mesh_{na}x{nc}"

    # One common policy for every mesh.  It brackets Nguyen et al.'s 80 equal
    # increments while allowing automatic growth/shrinkage.
    settings = SolverSettings(
        initial_step=.0125, maximum_step=.025, minimum_step=.000390625,
        residual_tolerance=1e-4, maximum_newton_iterations=32,
        eigenvalue_relative_floor=1e-11)
    diagnostic_loads = (.25, .5, .75, 1.)

    def accepted(load, values, result, step):
        a, b, c = problem.observables(values)
        parts = problem.parts(values)
        row = {
            "load_factor": float(load), "full_force_each": float(load*PMAX),
            "reduced_applied_force": float(load*REDUCED_FORCE),
            "w_A": a, "minus_u_B": b, "minus_u_C": c,
            "newton_iterations": result["iterations"],
            "minimum_line_search_alpha": result["minimum_alpha"],
            "accepted_step_size": float(step), "solver_status": result["status"],
            "residual_inf": result["history"][-1]["residual_inf"],
        }
        row.update(_energy_fields(parts))
        return row

    continuation = adaptive_continuation(
        lambda load, _anchor: lambda x: problem.potential(x, load),
        problem.initial, accepted, settings, mandatory_loads=diagnostic_loads)
    rows = continuation["accepted"]
    completed = continuation["status"] == "CONVERGED"
    metrics = _path_metrics(rows) if completed else None
    meaningful = [row for row in rows if abs(row["internal_energy"]) > 1e-10]
    max_xm = max((row["XM_fraction_internal"] for row in meaningful), default=None)
    max_xw = max((row["XW_fraction_internal"] for row in meaningful), default=None)
    summary = {
        "case": "open_cylinder_pullout",
        "classification": (
            "OPEN_CYLINDER_PULLOUT_VALIDATED_WITH_FINITE_MESH_ERROR"
            if completed else
            "OPEN_CYLINDER_PULLOUT_SOLVER_COMPLEXITY_NOT_WORTH_DELAY"),
        "mesh_cells_axial": na, "mesh_cells_circumferential": nc,
        "vertices": (na+1)*(nc+1), "scalar_dof": 3*(na+1)*(nc+1),
        "active_dof": problem.initial.numel(),
        "parameters": {"L": LENGTH, "R": RADIUS, "t": THICKNESS,
                       "E": E, "nu": NU, "Pmax_each": PMAX,
                       "reduced_force": REDUCED_FORCE},
        "domain": "one eighth: x in [0,L/2], theta in [0,pi/2]",
        "boundary_condition": (
            "three physical-midpoint reflection edges; x=L/2 natural free edge; "
            "raw gauge quotient slice separate"),
        "load_potential": "lambda*(Pmax/4)*Delta z_A; global dead load",
        "observables": "w_A=Delta z_A, -u_B=-Delta y_B, -u_C=-Delta y_C",
        "formulation_fingerprint": FORMULATION_FINGERPRINT,
        "reference_status": "DIGITIZED_REFERENCE",
        "reference_note": "Sze markers and Duong quartic vector curve, Fig. 12c",
        "reduction_audit": {key: value for key, value in problem.reduction.items()
                            if key not in ("transform", "point_a", "point_b", "point_c")},
        "small_load_gate": gate, "continuation_status": continuation["status"],
        "last_converged_load": continuation["last_converged_load"],
        "arc_length_used": False, "path_metrics": metrics,
        "maximum_meaningful_XM_fraction_internal": max_xm,
        "maximum_meaningful_XW_fraction_internal": max_xw,
        "solver_settings": continuation["settings"],
        "steps": rows, "rejected_steps": continuation["rejected"],
    }
    write_json(destination/f"{stem}.json", summary)
    write_csv(destination/f"{stem}_raw_path.csv", rows)
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--meshes", nargs="+", default=("8x8", "12x12", "16x16"))
    parser.add_argument("--skip-small-load-gate", action="store_true")
    args = parser.parse_args()
    results = []
    for mesh in args.meshes:
        na, nc = (int(value) for value in mesh.lower().split("x"))
        result = solve(na, nc, run_small_gate=not args.skip_small_load_gate)
        results.append(result)
        if result["classification"].endswith("NOT_WORTH_DELAY"):
            break
    write_json(CASE_OUTPUT/"results.json", results)
    write_csv(CASE_OUTPUT/"summary.csv", [{
        "mesh": f'{row["mesh_cells_axial"]}x{row["mesh_cells_circumferential"]}',
        "classification": row["classification"],
        "last_converged_load": row.get("last_converged_load"),
        "small_load_gate": row.get("small_load_gate", {}).get("passed"),
        "path_RMS": ((row.get("path_metrics") or {}).get("combined") or {}).get(
            "RMS_normalized_path_error"),
    } for row in results])
    write_json(CASE_OUTPUT/"summary.json", results)
    print(results[-1]["classification"], results[-1].get("last_converged_load"))


if __name__ == "__main__":
    main()
