"""End-moment cantilever benchmark on the current production energy.

Ported from ``CB-shell/experiments/checkerboard_nonlinear_benchmarks/
cantilever.py`` at tag ``v1.0-arxiv``.  The geometry, load potential, and
kinematic reduction are unchanged; the energy/targets are the current
``checkerboard_shell`` production functions.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import torch

from checkerboard_shell.mechanics.energy import make_targets
from experiments.cmame_round4.common import energy_parts, write_csv, write_json
from experiments.cmame_round4.nonlinear_solver import (
    SolverSettings, adaptive_continuation, interpolate_curve, solve_load_level,
)


YOUNG = 1.2e6
POISSON = 0.0
LENGTH = 12.0
WIDTH = 1.0
THICKNESS = 0.1
MAXIMUM_MOMENT = 50.0*math.pi/3.0


def analytical_tip(load_factor: float) -> tuple[float, float, float]:
    angle = 2.0*math.pi*float(load_factor)
    if abs(angle) < 1e-12:
        return LENGTH, 0.0, angle
    return (LENGTH*math.sin(angle)/angle,
            LENGTH*(1.0-math.cos(angle))/angle, angle)


def reference_grid(nx: int, ny: int, dtype=torch.float64):
    x = torch.linspace(0.0, LENGTH, nx+1, dtype=dtype)
    y = torch.linspace(-WIDTH/2.0, WIDTH/2.0, ny+1, dtype=dtype)
    yy, xx = torch.meshgrid(y, x, indexing="ij")
    q = torch.stack((xx, yy, torch.zeros_like(xx)), -1)
    return q, LENGTH/nx, WIDTH/ny


class CantileverProblem:
    def __init__(self, nx: int, ny: int = 4):
        if ny < 3:
            raise ValueError("at least three cells across the width are required")
        self.nx, self.ny = nx, ny
        self.reference, self.hx, self.hy = reference_grid(nx, ny)
        self.targets = make_targets(self.reference, self.hx, self.hy)
        # The published pure-bending solution is the plane-deformation branch
        # (nu=0).  Four cells are retained across the width so every frozen
        # B/W boundary stencil exists, while one x/z pair is shared by all
        # width rows.  This is a kinematic symmetry reduction, not a penalty.
        self.initial = torch.zeros((nx-1)*2, dtype=self.reference.dtype)

    def current(self, values: torch.Tensor) -> torch.Tensor:
        free_displacement = values.reshape(self.nx-1, 2)
        current = self.reference.clone()
        current[:, 2:, 0] += free_displacement[:, 0]
        current[:, 2:, 2] += free_displacement[:, 1]
        # Exact one-sided boundary position and tangent constraints on each
        # width line: -3 u_0+4 u_1-u_2=0 with u_0=0.
        current[:, 1, 0] += 0.25*free_displacement[0, 0]
        current[:, 1, 2] += 0.25*free_displacement[0, 1]
        return current

    @staticmethod
    def physical_centerline(vertices: torch.Tensor) -> torch.Tensor:
        """Gauge-invariant longitudinal centerline from width-edge midpoints."""
        width_midpoints = 0.5*(vertices[:-1] + vertices[1:])
        return width_midpoints.mean(0)

    @staticmethod
    def turning_angle(vertices: torch.Tensor) -> torch.Tensor:
        center = CantileverProblem.physical_centerline(vertices)
        tangents = center[1:] - center[:-1]
        cross = (tangents[:-1, 0]*tangents[1:, 2]
                 - tangents[:-1, 2]*tangents[1:, 0])
        dot = (tangents[:-1, 0]*tangents[1:, 0]
               + tangents[:-1, 2]*tangents[1:, 2])
        return torch.atan2(cross, dot).sum()

    @staticmethod
    def end_tangent(vertices: torch.Tensor) -> torch.Tensor:
        center = CantileverProblem.physical_centerline(vertices)
        return 3.0*center[-1]-4.0*center[-2]+center[-3]

    def moment_work(self, values: torch.Tensor, load_factor: float,
                    anchor_tangent: torch.Tensor,
                    anchor_work: torch.Tensor | float = 0.0) -> torch.Tensor:
        tangent = self.end_tangent(self.current(values))
        anchor = anchor_tangent.to(tangent)
        cross = anchor[0]*tangent[2]-anchor[2]*tangent[0]
        dot = anchor[0]*tangent[0]+anchor[2]*tangent[2]
        local_rotation = torch.atan2(cross, dot)
        return float(load_factor)*MAXIMUM_MOMENT*local_rotation + anchor_work

    def parts(self, values: torch.Tensor) -> dict:
        current = self.current(values)
        parts = energy_parts(
            current, self.reference, self.hx, self.hy, self.targets,
            YOUNG, POISSON, THICKNESS)
        parts["_vertices"] = current
        return parts

    def potential(self, values: torch.Tensor, load_factor: float,
                  anchor_tangent: torch.Tensor | None = None,
                  anchor_work: torch.Tensor | float = 0.0):
        current = self.current(values)
        parts = energy_parts(
            current, self.reference, self.hx, self.hy, self.targets,
            YOUNG, POISSON, THICKNESS)
        if anchor_tangent is None:
            anchor_tangent = self.end_tangent(self.reference)
        work = self.moment_work(values, load_factor, anchor_tangent, anchor_work)
        return parts["total"]-work

    def tip(self, values: torch.Tensor):
        tip = self.physical_centerline(self.current(values))[-1]
        return float(tip[0]), float(tip[2])


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


def solve(nx: int, ny: int = 4, output: Path | None = None,
          run_small_gate: bool = False) -> dict:
    problem = CantileverProblem(nx, ny)
    gate = {"passed": None, "rows": []}
    del run_small_gate  # the archived small-load gate is not rerun here

    settings = SolverSettings(initial_step=0.025, maximum_step=0.05,
                              minimum_step=0.000390625,
                              residual_tolerance=5e-7,
                              maximum_newton_iterations=45)

    work_offset = 0.0

    def make_potential(load, accepted):
        nonlocal work_offset
        anchor = problem.end_tangent(problem.current(accepted)).detach()
        previous_angle = problem.turning_angle(
            problem.current(accepted)).detach()
        work_offset = float(load)*MAXIMUM_MOMENT*previous_angle
        return lambda values: problem.potential(
            values, load, anchor, work_offset)

    def on_accept(load, values, result, step):
        parts = problem.parts(values)
        tip_x, tip_z = problem.tip(values)
        exact_x, exact_z, exact_angle = analytical_tip(load)
        row = {
            "load_factor": float(load),
            "moment": float(load*MAXIMUM_MOMENT),
            "tip_x": tip_x, "tip_z": tip_z,
            "tip_u": tip_x-LENGTH, "tip_w": tip_z,
            "analytical_tip_x": exact_x, "analytical_tip_z": exact_z,
            "tip_vector_error": float(math.hypot(tip_x-exact_x, tip_z-exact_z)),
            "turning_angle": float(problem.turning_angle(
                problem.current(values))),
            "analytical_turning_angle": exact_angle,
            "newton_iterations": result["iterations"],
            "minimum_line_search_alpha": result["minimum_alpha"],
            "accepted_step_size": float(step),
            "solver_status": result["status"],
            "residual_inf": result["history"][-1]["residual_inf"],
        }
        row.update(_energy_fields(parts))
        return row

    continuation = adaptive_continuation(
        make_potential, problem.initial, on_accept, settings,
        mandatory_loads=(0.25, 0.5, 0.75, 1.0))
    rows = continuation["accepted"]
    sample = np.linspace(0.05, 1.0, 20)
    interpolated = interpolate_curve(
        rows, sample, ("tip_x", "tip_z", "turning_angle"))
    exact = np.asarray([analytical_tip(load) for load in sample])
    numerical = np.stack((interpolated["tip_x"], interpolated["tip_z"]), -1)
    difference = numerical-exact[:, :2]
    scale = np.linalg.norm(exact[:, :2])
    summary = {
        "case": "cantilever_under_end_moment",
        "classification": ("VALIDATED"
                           if continuation["status"] == "CONVERGED"
                           else "SOLVER_LIMITED"),
        "mesh_cells": [nx, ny],
        "vertices": (nx+1)*(ny+1),
        "scalar_dof": 3*(nx+1)*(ny+1),
        "active_dof": problem.initial.numel(),
        "parameters": {"E": YOUNG, "nu": POISSON, "L": LENGTH,
                       "W": WIDTH, "t": THICKNESS,
                       "Mmax": MAXIMUM_MOMENT},
        "formulation_fingerprint": {
            "first_form": "a_B=T_B^T T_B",
            "membrane_energy": "t/4*[sum_B A_B^0 Q(Delta a_B)+8*mu*X_M^2]",
            "second_form": "accepted symmetric b_W",
            "bending_compatibility": "accepted objective/reference-relative X_W_rel",
            "state": "positions only modulo exact checkerboard gauge",
            "energy_source": "checkerboard_shell.mechanics.energy.paper_components",
        },
        "boundary_observation": (
            "gauge-invariant centerline of physical width-edge midpoints"),
        "symmetry_reduction": (
            "published nu=0 plane-deformation branch; four width cells retain "
            "complete frozen B/W stencils"),
        "continuation_status": continuation["status"],
        "small_load_gate": gate,
        "maximum_turning_angle": rows[-1]["turning_angle"] if rows else 0.0,
        "path_RMS_normalized_error": float(np.linalg.norm(difference)/scale),
        "path_max_normalized_error": float(
            np.linalg.norm(difference, axis=1).max()/LENGTH),
        "final_tip_error": (rows[-1]["tip_vector_error"]/LENGTH
                            if rows else None),
        "steps": rows,
        "rejected_steps": continuation["rejected"],
        "arc_length_used": False,
    }
    destination = output or Path(__file__).resolve().parent/"raw_data"/"cantilever"
    write_json(destination/f"mesh_{nx}x{ny}.json", summary)
    write_csv(destination/f"mesh_{nx}x{ny}_raw_steps.csv", rows)
    return summary


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--meshes", nargs="+", type=int,
                        default=(16, 24, 32))
    parser.add_argument("--width-cells", type=int, default=4)
    arguments = parser.parse_args()
    payload = [solve(nx, arguments.width_cells) for nx in arguments.meshes]
    print([(p["mesh_cells"], p["classification"],
            p.get("path_RMS_normalized_error")) for p in payload])
