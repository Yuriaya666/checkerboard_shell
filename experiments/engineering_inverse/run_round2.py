"""One-command rebuild of the mechanics-consistency follow-up.

Full execution is CPU intensive; each stage writes a resumable result file.
Use ``--quick`` for tests, fixed-input diagnostics, saved-design references
and figures without rerunning the long nonlinear optimizations/search.
"""
from __future__ import annotations

import argparse
import json
import time

from experiments.engineering_inverse import (
    aircraft_nonlinear_reference, aircraft_meshcheck, automotive_spring_design,
    automotive_spring_reference, automotive_mitc4_spring_design,
    automotive_corotational_design, beam_kinematics_audit,
    diagnose_reactions, mount_springs, nonlinear_audit,
    performance_round2, make_round2_figures, gradient_cost_round2,
    automotive_design_crosscheck, profile_round2)
from experiments.engineering_inverse.aircraft_skin import OUT


def run(full=True):
    stages=[]
    def stage(name,call):
        started=time.perf_counter()
        call()
        stages.append(dict(name=name,seconds=time.perf_counter()-started))
        (OUT/"round2_execution_log.json").write_text(json.dumps(stages,indent=2)+"\n")
    stage("hard_patch_fixed_input_diagnosis",diagnose_reactions.run)
    stage("area_integrated_bracket_fixed_input_audit",mount_springs.audit)
    if full:
        stage("checkerboard_nonlinear_bracket_design",automotive_spring_design.nonlinear_design)
    stage("automotive_common_reference",automotive_spring_reference.run)
    if full:
        stage("mitc4_linear_bracket_design",automotive_mitc4_spring_design.run)
        stage("corotational_bracket_refinement",automotive_corotational_design.run)
        stage("same_model_gradient_cost",gradient_cost_round2.run)
        stage("both_automotive_designs_common_reference",automotive_design_crosscheck.run)
    stage("aircraft_beam_kinematics",beam_kinematics_audit.run)
    if full:
        stage("aircraft_low_cardinality_nonlinear_audit",nonlinear_audit.run)
    stage("aircraft_checkerboard_meshcheck",aircraft_meshcheck.run)
    stage("aircraft_corotational_reference",aircraft_nonlinear_reference.run)
    if full:
        stage("aircraft_31_layout_performance",performance_round2.measure)
        stage("aircraft_operator_profile",profile_round2.run)
    stage("paper_vector_figures",lambda:(make_round2_figures.automotive(),
                                          make_round2_figures.aircraft()))
    return stages


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--quick",action="store_true")
    args=parser.parse_args()
    print(json.dumps(run(full=not args.quick),indent=2))
