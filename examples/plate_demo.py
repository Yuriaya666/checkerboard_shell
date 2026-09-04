"""Minimal rectangular Checkerboard Shell forward solve.

Run from the project directory with ``python examples/plate_demo.py``.
"""

from pathlib import Path
import sys

# Keep the source checkout directly runnable before an editable install.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import torch

from checkerboard_shell import (
    BoundaryConditions, CheckerboardMesh, CheckerboardShell,
    IsotropicMaterial, NodalLoad, NewtonSettings,
)


def main():
    mesh = CheckerboardMesh.rectangle(width=1.0, height=0.6, nx=6, ny=4)
    material = IsotropicMaterial(young=2.0e5, poisson=0.0, thickness=0.04)
    boundary = BoundaryConditions.clamped_boundary(mesh, "left")
    tip_nodes = mesh.boundaries["right"]
    load = NodalLoad.total_force(
        mesh, force=(0.0, 0.0, -1.0e-1), node_ids=tip_nodes,
        name="downward_tip_resultant")
    settings = NewtonSettings(
        residual_tolerance=1.0e-9, maximum_newton_iterations=20)
    model = CheckerboardShell(mesh, material, boundary, load, settings)
    result = model.solve()

    print(f"status: {result.status}")
    print(f"iterations: {result.iterations}")
    print(f"residual_inf: {result.residual_inf:.3e}")
    print(f"max displacement: {result.max_displacement:.6e} m")
    print("deformed vertices (first five):")
    print(result.deformed_vertices.reshape(-1, 3)[:5])
    # Plotting is optional; save a portable image when Matplotlib is present.
    result.plot(save_path=Path("outputs/plate_demo.png"), displacement_scale=1.0)
    result.save(Path("outputs/plate_demo.npz"))
    return result


if __name__ == "__main__":
    main()
