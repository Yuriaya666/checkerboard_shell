"""Run all four design demos and write a compact summary."""

from __future__ import annotations

import json

from experiments.cmame_round5 import (
    actuation_design, adjoint_check, large_rotation_design, shape_design,
)
from experiments.cmame_round5.common import DATA, LOGS, write_json


def run() -> dict:
    adjoint = adjoint_check.run()
    actuation = actuation_design.run()
    shape = shape_design.run()
    rotation = large_rotation_design.run()
    summary = {
        "adjoint_max_relative_error": adjoint["maximum_relative_error"],
        "actuation": {
            "displacement_relative_error": actuation["displacement_relative_error"],
            "force_relative_error": actuation["force_relative_error"],
            "iterations": actuation["iterations"],
            "design_variables": actuation["design_variables"],
            "free_dofs": actuation["free_dofs"],
        },
        "shape": {
            "displacement_relative_error": shape["displacement_relative_error"],
            "height_relative_error": shape["height_relative_error"],
            "iterations": shape["iterations"],
            "design_variables": shape["design_variables"],
            "free_dofs": shape["free_dofs"],
        },
        "large_rotation": {
            "tip_error": rotation["tip_error"],
            "iterations": rotation["iterations"],
            "design_variables": rotation["design_variables"],
            "free_dofs": rotation["free_dofs"],
        },
    }
    write_json(DATA/"design_summary.json", summary)
    write_json(LOGS/"run_all_summary.json", summary)
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
