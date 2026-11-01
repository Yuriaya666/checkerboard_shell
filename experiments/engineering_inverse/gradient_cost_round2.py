"""Measured same-model implicit versus centered-FD automotive design gradient."""
from __future__ import annotations

import json
import time
import numpy as np
import torch

from experiments.engineering_inverse.automotive_panel import build_case
from experiments.engineering_inverse.mount_springs import checkerboard_forward,bracket_matrices
from experiments.engineering_inverse.automotive_spring_design import response_sensitivity
from experiments.engineering_inverse.aircraft_skin import OUT


def run(step_m=1e-5):
    saved=json.loads((OUT/"automotive_spring_design_12x8.json").read_text())
    x=np.array(saved["adjustment_mm"])
    d=x*1e-3;case=build_case()
    weight=np.tile([8.,1.],case["mesh"].nx+1)
    start=time.perf_counter()
    base=checkerboard_forward(case,d)
    base_seconds=time.perf_counter()-start
    if base["status"]!="CONVERGED":raise RuntimeError("baseline did not converge")
    error=(case["geometry_error"]+case["O"]@base["displacement"].ravel())*1e3
    K,D,W=bracket_matrices(case)
    start=time.perf_counter()
    J,dF=response_sensitivity(case,base,K,D,W)
    implicit_seconds=time.perf_counter()-start
    exact=2*J.T@(weight*error)/len(error)+.001*x
    initial=torch.as_tensor(base["displacement"]).reshape(case["mesh"].shape)
    def objective(input_m):
        s=checkerboard_forward(case,input_m,initial=initial)
        if s["status"]!="CONVERGED":raise RuntimeError(s["status"])
        e=(case["geometry_error"]+case["O"]@s["displacement"].ravel())*1e3
        y=input_m*1e3
        return float((weight*e)@e/len(e)+.0005*np.dot(y,y))
    start=time.perf_counter()
    fd=np.empty(8)
    for k in range(8):
        perturb=np.zeros(8);perturb[k]=step_m
        fd[k]=(objective(d+perturb)-objective(d-perturb))/(2*step_m*1e3)
    fd_seconds=time.perf_counter()-start
    report=dict(mesh=[12,8],baseline_seconds=base_seconds,
        implicit_linearization_seconds=implicit_seconds,
        finite_difference_seconds=fd_seconds,
        finite_difference_forward_equilibria=16,
        step_m=step_m,implicit_gradient=exact.tolist(),
        finite_difference_gradient=fd.tolist(),
        max_absolute_gradient_difference=float(np.max(np.abs(exact-fd))),
        max_relative_to_gradient_scale=float(np.max(np.abs(exact-fd))/
            max(np.max(abs(exact)),np.max(abs(fd)),1e-12)))
    (OUT/"gradient_cost_round2.json").write_text(json.dumps(report,indent=2)+"\n")
    return report


if __name__=="__main__":print(run())
