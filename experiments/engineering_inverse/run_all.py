"""One-command reproducible run: python -m experiments.engineering_inverse.run_all.

The full run includes all 31 linear layouts and 15 nonlinear low-cardinality
audits. It can take tens of minutes on a single CPU core.
"""

import argparse
import json
import platform
import time

import numpy as np
import scipy
import torch

from experiments.engineering_inverse import (
    aircraft_skin,automotive_panel,automotive_reference,
    automotive_crosscheck,reference_check,make_figures,nonlinear_audit,
    three_workflows,timing_study,load_step_check)


def run(full=True):
    log=dict(started_utc=__import__("datetime").datetime.now(
        __import__("datetime").timezone.utc).isoformat(),
        environment=dict(python=platform.python_version(),numpy=np.__version__,
                         scipy=scipy.__version__,torch=torch.__version__,
                         processor=platform.processor()),stages=[])
    def stage(name,fn):
        start=time.perf_counter()
        fn()
        log["stages"].append(dict(name=name,seconds=time.perf_counter()-start))
        (aircraft_skin.OUT/"run_log.json").write_text(json.dumps(log,indent=2)+"\n")
    stage("aircraft_linear_search",lambda:aircraft_skin.run(
        *aircraft_skin.PARAMS["primary_mesh"]))
    stage("automotive_linear_and_nonlinear",lambda:automotive_panel.run(
        *automotive_panel.PARAMS["primary_mesh"]))
    if full:
        stage("aircraft_mesh_screens",lambda:(aircraft_skin.run(16,32,True),
                                              aircraft_skin.run(20,40,True)))
        stage("automotive_mesh_screens",lambda:(automotive_panel.run(16,12,True,False),
                                                automotive_panel.run(20,16,True,False)))
        stage("aircraft_low_cardinality_nonlinear_audit",nonlinear_audit.run)
        def final_design():
            result=json.loads((aircraft_skin.OUT/"aircraft_12x24.json").read_text())
            case=aircraft_skin.build_case()
            selected=aircraft_skin.PARAMS["nonlinear_verified_layout_indices"]
            row=next(r for r in result["checkerboard"]["search"]
                     if r["stations"]==selected)
            modes=[aircraft_skin.nonlinear_optimize_controls(
                case,selected,row["cases"][k]["strokes_m"],case["targets"][k],
                max_iterations=4) for k in range(2)]
            payload=dict(stations=selected,cases=modes,
                         all_feasible=all(m["success"] for m in modes),
                         source="exact-tangent nonlinear coupled shell + position beams")
            (aircraft_skin.OUT/"aircraft_nonlinear_final_12x24.json").write_text(
                json.dumps(payload,indent=2)+"\n")
        stage("aircraft_final_nonlinear_design",final_design)
        stage("aircraft_load_step_check",load_step_check.run)
        stage("aircraft_common_reference",reference_check.run)
        stage("automotive_common_reference",automotive_reference.run)
        stage("automotive_crosscheck",automotive_crosscheck.run)
        stage("three_aircraft_workflows",three_workflows.run)
        stage("timing_repeats",timing_study.run)
        stage("vector_figures_and_csv",lambda:(make_figures.aircraft(),
                                               make_figures.automotive()))
    return log


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--quick",action="store_true",
                        help="run only primary linear screens and automotive nonlinear check")
    args=parser.parse_args()
    print(json.dumps(run(not args.quick),indent=2))
