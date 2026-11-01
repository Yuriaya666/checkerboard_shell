"""Aircraft case-2 load/stroke continuation against saved direct equilibrium."""

import json
import time

import numpy as np
import scipy.sparse as sp
import torch

from checkerboard_shell import BoundaryConditions,CheckerboardShell,NewtonSettings
from experiments.engineering_inverse.aircraft_skin import OUT,build_case
from experiments.engineering_inverse.common import positional_beam_stiffness


def run():
    case=build_case()
    final=json.loads((OUT/"aircraft_nonlinear_final_12x24.json").read_text())
    selected=final["stations"]
    d=np.asarray(final["cases"][1]["strokes_m"])
    groups=[[3*n+2 for n in case["patches"][k]] for k in selected]
    fixed=sorted(set(case["spar"]+[dof for group in groups for dof in group]))
    mask=torch.zeros(case["mesh"].shape,dtype=torch.bool)
    mask.reshape(-1)[fixed]=True
    Kbeam=sum((positional_beam_stiffness(case["mesh"].reference_vertices.numpy(),*line)
               for line in case["beam_lines"]),
              np.zeros((3*case["mesh"].num_vertices,)*2))
    prescribed=torch.zeros(case["mesh"].shape,dtype=torch.float64)
    model=CheckerboardShell(case["mesh"],case["material"],
        BoundaryConditions(mask,prescribed),case["load"],
        NewtonSettings(residual_tolerance=1e-5,maximum_newton_iterations=30),
        sp.csr_matrix(Kbeam))
    state=None;rows=[]
    start=time.perf_counter()
    for fraction in (.25,.5,.75,1.0):
        values=torch.zeros(case["mesh"].shape,dtype=torch.float64)
        for group,stroke in zip(groups,d):values.reshape(-1)[group]=fraction*stroke
        model.boundary_conditions=BoundaryConditions(mask,values)
        solved=model.solve(initial_displacement=state,load_factor=fraction)
        state=(solved.deformed_vertices-case["mesh"].reference_vertices).detach()
        rows.append(dict(fraction=fraction,status=solved.status,
                         newton_iterations=solved.iterations,
                         residual_inf=solved.residual_inf))
        if solved.status!="CONVERGED":break
    reference=np.asarray(final["cases"][1]["displacement_m"])
    difference=state.numpy().reshape(-1)-reference
    data=dict(rows=rows,all_converged=all(r["status"]=="CONVERGED" for r in rows),
        maximum_position_difference_mm=float(np.max(np.abs(difference))*1e3),
        rms_position_difference_mm=float(np.sqrt(np.mean(difference**2))*1e3),
        elapsed_seconds=time.perf_counter()-start)
    (OUT/"aircraft_load_step_check.json").write_text(json.dumps(data,indent=2)+"\n")
    return data


if __name__=="__main__":
    print(json.dumps(run(),indent=2))
