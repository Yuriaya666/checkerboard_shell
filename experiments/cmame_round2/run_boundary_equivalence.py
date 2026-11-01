"""Old dense versus new local-stencil energy/derivative equivalence."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import torch

from checkerboard_shell.material.isotropic import IsotropicMaterial
from checkerboard_shell.mechanics.energy import make_targets, paper_components


ROOT=Path(__file__).resolve().parent
OUT=ROOT/"01_sparse_boundary"/"raw_data"
HIST=Path("/home/song/projects/CB-shell")
OLD=OUT/"old_dense_derivatives.npz"

PROGRAM=r'''
import numpy as np, torch, sys
from experiments.checkerboard_benchmark.mechanics.checkerboard_koiter import koiter_targets
from experiments.checkerboard_paper_benchmark.formulation import paper_components
n=4;theta0=np.deg2rad(18.)
theta=torch.linspace(theta0,np.pi/2,n+1,dtype=torch.float64);phi=torch.linspace(0,np.pi/2,n+1,dtype=torch.float64)
tt,pp=torch.meshgrid(theta,phi,indexing='ij');ref=10*torch.stack((torch.sin(tt)*torch.cos(pp),torch.sin(tt)*torch.sin(pp),torch.cos(tt)),-1)
hx=(np.pi/2)/n;hy=(np.pi/2-theta0)/n;targets=koiter_targets(ref,hx,hy)
phase=torch.arange(ref.numel(),dtype=ref.dtype).reshape_as(ref);q=2e-3*torch.sin(.29*phase);flat=q.reshape(-1).requires_grad_(True)
def energy(x):return paper_components(ref+x.reshape_as(ref),ref,hx,hy,targets,6.825e7,.3,.04)['total']
e=energy(flat);r,=torch.autograd.grad(e,flat);k=torch.autograd.functional.hessian(energy,flat,vectorize=True)
np.savez(sys.argv[1],energy=float(e.detach()),residual=r.detach().numpy(),tangent=k.detach().numpy())
'''


def run():
    OUT.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy();env["PYTHONPATH"]=str(HIST)
    subprocess.run([sys.executable,"-c",PROGRAM,str(OLD)],cwd=HIST,env=env,check=True)
    old=np.load(OLD);n=4;theta0=np.deg2rad(18.)
    theta=torch.linspace(theta0,np.pi/2,n+1,dtype=torch.float64);phi=torch.linspace(0,np.pi/2,n+1,dtype=torch.float64)
    tt,pp=torch.meshgrid(theta,phi,indexing="ij");ref=10*torch.stack((torch.sin(tt)*torch.cos(pp),torch.sin(tt)*torch.sin(pp),torch.cos(tt)),-1)
    hx=(np.pi/2)/n;hy=(np.pi/2-theta0)/n;targets=make_targets(ref,hx,hy);material=IsotropicMaterial(6.825e7,.3,.04)
    phase=torch.arange(ref.numel(),dtype=ref.dtype).reshape_as(ref);q=2e-3*torch.sin(.29*phase);flat=q.reshape(-1).requires_grad_(True)
    def energy(x):return paper_components(ref+x.reshape_as(ref),ref,hx,hy,targets,material)["total"]
    e=energy(flat);r,=torch.autograd.grad(e,flat);k=torch.autograd.functional.hessian(energy,flat,vectorize=True)
    new_e=float(e.detach());new_r=r.detach().numpy();new_k=k.detach().numpy()
    payload={"case":"finite-deformed curved M3 n=4","energy_absolute_difference":abs(new_e-float(old["energy"])),
        "residual_relative_difference":float(np.linalg.norm(new_r-old["residual"])/np.linalg.norm(old["residual"])),
        "tangent_relative_difference":float(np.linalg.norm(new_k-old["tangent"])/np.linalg.norm(old["tangent"])),
        "new_energy":new_e,"old_energy":float(old["energy"])}
    (OUT/"boundary_equivalence.json").write_text(json.dumps(payload,indent=2)+"\n")
    return payload


if __name__=="__main__":print(json.dumps(run(),indent=2))

