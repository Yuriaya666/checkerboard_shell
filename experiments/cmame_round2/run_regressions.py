"""Post-optimization consistency and gauge-invariance regression suite."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import torch

from checkerboard_shell.geometry.mesh import CheckerboardMesh
from checkerboard_shell.material.isotropic import IsotropicMaterial
from checkerboard_shell.mechanics.energy import make_targets, paper_components
from checkerboard_shell.mechanics.low_rank import assemble_paper_low_rank_tangent
from experiments.cmame_additional.common import invariant_metrics
from experiments.cmame_round2.recovered_drivers import m3_problem, solve_curved_sparse


ROOT = Path(__file__).resolve().parent


def json_finite(value):
    """Represent undefined relative errors as JSON null, not nonstandard NaN."""
    if isinstance(value, dict):
        return {key: json_finite(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_finite(item) for item in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def run() -> dict:
    direction = torch.tensor([1.0, 0.7, -0.4], dtype=torch.float64)
    direction /= torch.linalg.vector_norm(direction)
    multipliers = (0.0, 1e-8, 1e-4, 1e-2, 0.1, 0.5, 1.0)
    states = []
    m3 = m3_problem(6)
    states.append(("undeformed_curved_M3", m3["reference"], m3["reference"],
                   m3["hx"], m3["hy"], m3["material"]))
    displacement, _ = solve_curved_sparse(m3)
    current_m3 = m3["reference"]+torch.as_tensor(
        displacement, dtype=m3["reference"].dtype).reshape_as(m3["reference"])
    states.append(("linear_M3_deformed", m3["reference"], current_m3,
                   m3["hx"], m3["hy"], m3["material"]))
    archive = np.load(ROOT.parent/"cmame_additional"/"02_gauge_invariance"/
                      "raw_data"/"state_nonlinear_cantilever.npz")
    reference = torch.as_tensor(archive["reference"], dtype=torch.float64)
    current = torch.as_tensor(archive["current"], dtype=torch.float64)
    states.append(("nonlinear_cantilever", reference, current, 2.0/6, 0.6/4,
                   IsotropicMaterial(2e4, 0.3, 0.08)))
    gauge_rows = []
    for name, reference, current, hx, hy, material in states:
        length = float(torch.linalg.vector_norm(
            reference.reshape(-1,3).max(0).values-reference.reshape(-1,3).min(0).values))
        rows = invariant_metrics(reference,current,hx,hy,material,
                                 [a*length for a in multipliers],direction)
        for multiplier,row in zip(multipliers,rows):
            gauge_rows.append({"state":name,"alpha_over_L":multiplier,**row})

    # Independent objectivity/tangent finite-difference checks on a deformed
    # curved patch, using the optimized local stencil and low-rank path.
    reference=m3["reference"]; material=m3["material"]; targets=m3["targets"]
    phase=torch.arange(reference.numel(),dtype=reference.dtype).reshape_as(reference)
    q=1e-3*torch.sin(.19*phase)
    parts=paper_components(reference+q,reference,m3["hx"],m3["hy"],targets,material)
    translation=torch.tensor([.31,-.22,.17],dtype=reference.dtype)
    angle=.37; rotation=torch.tensor([[np.cos(angle),-np.sin(angle),0],
        [np.sin(angle),np.cos(angle),0],[0,0,1]],dtype=reference.dtype)
    translated=paper_components(reference+q+translation,reference+translation,
        m3["hx"],m3["hy"],make_targets(reference+translation,m3["hx"],m3["hy"]),material)
    rotated_ref=reference@rotation.T; rotated=(reference+q)@rotation.T
    rotated_parts=paper_components(rotated,rotated_ref,m3["hx"],m3["hy"],
        make_targets(rotated_ref,m3["hx"],m3["hy"]),material)
    tangent,residual,_=assemble_paper_low_rank_tangent(q,reference,m3["hx"],m3["hy"],targets,material)
    vector=np.cos(np.arange(reference.numel())*.41); eps=1e-7
    def residual_at(flat):
        value=flat.detach().clone().requires_grad_(True)
        energy=paper_components(reference+value.reshape_as(reference),reference,
            m3["hx"],m3["hy"],targets,material)["total"]
        grad,=torch.autograd.grad(energy,value);return grad.detach().numpy()
    fd=(residual_at((q.reshape(-1)+eps*torch.as_tensor(vector)))
        -residual_at((q.reshape(-1)-eps*torch.as_tensor(vector))))/(2*eps)
    base_energy = abs(float(parts["total"]))
    translation_difference = abs(float(parts["total"]-translated["total"]))
    rotation_difference = abs(float(parts["total"]-rotated_parts["total"]))
    fd_error = float(np.linalg.norm(tangent.action(vector)-fd)/np.linalg.norm(fd))
    checks={
        "reference_energy":float(paper_components(reference,reference,m3["hx"],m3["hy"],targets,material)["total"]),
        "reference_residual_inf":float(np.max(np.abs(assemble_paper_low_rank_tangent(
            torch.zeros_like(reference),reference,m3["hx"],m3["hy"],targets,material)[1]))),
        "objectivity_base_energy": base_energy,
        "translation_energy_difference": translation_difference,
        "translation_energy_relative_difference": translation_difference/max(base_energy, 1e-300),
        "rotation_energy_difference": rotation_difference,
        "rotation_energy_relative_difference": rotation_difference/max(base_energy, 1e-300),
        "residual_directional_derivative_relative_error": fd_error,
        "tangent_action_fd_relative_error": fd_error,
        "local_tangent_symmetry_inf":float(np.max(np.abs((tangent.local-tangent.local.T).data))) if (tangent.local-tangent.local.T).nnz else 0.0,
        "gauge_max_relative_tangent_change":max(r["relative_tangent_change"] for r in gauge_rows),
        "gauge_max_physical_midpoint_change":max(r["physical_midpoint_change"] for r in gauge_rows),
    }
    payload=json_finite({"checks":checks,"gauge_rows":gauge_rows})
    out=ROOT/"reports";out.mkdir(parents=True,exist_ok=True)
    (out/"regression_results.json").write_text(json.dumps(payload,indent=2,allow_nan=False)+"\n")
    return payload


if __name__=="__main__": print(json.dumps(run(),indent=2,allow_nan=False))
