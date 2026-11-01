"""Surface-preserving distortion on the recovered manuscript problems."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

from checkerboard_shell.material.isotropic import IsotropicMaterial
from checkerboard_shell.mechanics.energy import make_targets
from checkerboard_shell.mechanics.low_rank import LowRankTangent, assemble_paper_low_rank_tangent
from experiments.cmame_round2.recovered_drivers import (
    M3_REFERENCE, PLATE_REFERENCE, ROOF_REFERENCE, m3_problem, roof_problem,
)


ROOT = Path(__file__).resolve().parent
OUT = ROOT/"04_distortion_revalidated"
SEED = 20260922
LEVELS = (0.0, 0.05, 0.10, 0.20, 0.30)
MESHES = {"plate": 16, "roof": 16, "M3": 12}


def logical_grid(n: int, delta: float, kind: str, benchmark: str):
    values = np.linspace(0.0, 1.0, n+1)
    eta, xi = np.meshgrid(values, values, indexing="ij")
    if delta == 0.0:
        return xi, eta
    if kind == "smooth":
        fx = np.sin(2*np.pi*xi)*np.sin(np.pi*eta)
        fy = np.sin(np.pi*xi)*np.sin(2*np.pi*eta)
    elif kind == "random":
        rng = np.random.default_rng(SEED+{"plate": 1, "roof": 2, "M3": 3}[benchmark])
        fx, fy = rng.normal(size=xi.shape), rng.normal(size=xi.shape)
        # Fixed, low-pass random tangential redistribution.
        for field in (fx, fy):
            for _ in range(2):
                padded = np.pad(field, 1, mode="edge")
                field[:] = sum(padded[a:a+n+1, b:b+n+1]
                               for a in range(3) for b in range(3))/9.0
        envelope = np.sin(np.pi*xi)*np.sin(np.pi*eta)
        fx *= envelope; fy *= envelope
        interior = np.s_[1:-1, 1:-1]
        target_x = np.sin(2*np.pi*xi)*np.sin(np.pi*eta)
        target_y = np.sin(np.pi*xi)*np.sin(2*np.pi*eta)
        fx *= np.sqrt(np.mean(target_x[interior]**2))/np.sqrt(np.mean(fx[interior]**2))
        fy *= np.sqrt(np.mean(target_y[interior]**2))/np.sqrt(np.mean(fy[interior]**2))
        if benchmark == "plate":
            fx[n//2, n//2] = 0.0; fy[n//2, n//2] = 0.0
    else:
        raise ValueError(kind)
    moved_x, moved_y = xi+(delta/n)*fx, eta+(delta/n)*fy
    moved_x[[0, -1], :] = xi[[0, -1], :]
    moved_x[:, [0, -1]] = xi[:, [0, -1]]
    moved_y[[0, -1], :] = eta[[0, -1], :]
    moved_y[:, [0, -1]] = eta[:, [0, -1]]
    return moved_x, moved_y


def param_areas(xi, eta):
    n = xi.shape[0]-1
    areas = np.zeros((n+1, n+1))
    signed = []
    for j in range(n):
        for i in range(n):
            x = np.array([xi[j,i], xi[j,i+1], xi[j+1,i+1], xi[j+1,i]])
            y = np.array([eta[j,i], eta[j,i+1], eta[j+1,i+1], eta[j+1,i]])
            area = 0.5*np.sum(x*np.roll(y,-1)-y*np.roll(x,-1))
            signed.append(area); areas[j:j+2, i:i+2] += area/4.0
    return areas, np.asarray(signed)


def quality(reference: torch.Tensor) -> dict:
    q = reference.numpy(); n = q.shape[0]-1
    scaled, aspect, angles = [], [], []
    for j in range(n):
        for i in range(n):
            p = [q[j,i], q[j,i+1], q[j+1,i+1], q[j+1,i]]
            normal = np.cross(p[1]-p[0], p[3]-p[0]); normal /= np.linalg.norm(normal)
            edges = [p[(k+1)%4]-p[k] for k in range(4)]
            lengths = [np.linalg.norm(v) for v in edges]; aspect.append(max(lengths)/min(lengths))
            for k in range(4):
                incoming = p[k]-p[(k-1)%4]; outgoing = p[(k+1)%4]-p[k]
                scaled.append(np.dot(np.cross(incoming,outgoing),normal)/(
                    np.linalg.norm(incoming)*np.linalg.norm(outgoing)))
                cosine = np.dot(-incoming,outgoing)/(np.linalg.norm(incoming)*np.linalg.norm(outgoing))
                angles.append(math.degrees(math.acos(np.clip(cosine,-1,1))))
    return {"minimum_scaled_jacobian": float(min(scaled)),
            "maximum_aspect_ratio": float(max(aspect)),
            "minimum_interior_angle_deg": float(min(angles)),
            "maximum_skew_angle_deg": float(max(abs(a-90) for a in angles))}


def reference_for(benchmark: str, n: int, xi, eta):
    u = torch.as_tensor(xi, dtype=torch.float64)
    v = torch.as_tensor(eta, dtype=torch.float64)
    if benchmark == "plate":
        return torch.stack((u, v, torch.zeros_like(u)), -1), 1/n, 1/n
    if benchmark == "roof":
        radius, half_angle = 25.0, math.radians(40.0)
        x = 25.0*u; theta = -half_angle+2*half_angle*v
        return torch.stack((x, radius*torch.sin(theta), radius*torch.cos(theta)), -1), 25/n, 2*radius*half_angle/n
    theta0 = math.radians(18.0)
    phi = 0.5*math.pi*u; theta = theta0+(0.5*math.pi-theta0)*v
    return 10*torch.stack((torch.sin(theta)*torch.cos(phi),
                           torch.sin(theta)*torch.sin(phi), torch.cos(theta)), -1), (0.5*math.pi)/n, (0.5*math.pi-theta0)/n


def solve_case(benchmark: str, n: int, xi, eta) -> tuple[dict, np.ndarray]:
    reference, hx, hy = reference_for(benchmark, n, xi, eta)
    targets = make_targets(reference, hx, hy)
    if benchmark == "plate":
        material = IsotropicMaterial(1.0, 0.3, 0.01)
    elif benchmark == "roof":
        material = IsotropicMaterial(4.32e8, 0.0, 0.25)
    else:
        material = IsotropicMaterial(6.825e7, 0.3, 0.04)
    tangent, residual, audit = assemble_paper_low_rank_tangent(
        torch.zeros_like(reference), reference, hx, hy, targets, material)
    weights, signed = param_areas(xi, eta)
    if np.min(signed) <= 0:
        raise ValueError("inverted logical cell")
    if benchmark == "plate":
        ids = []
        for j in range(1,n):
            for i in range(1,n): ids.append(3*(j*(n+1)+i)+2)
        active = np.asarray(ids); force = 1e-8*weights.reshape(-1)[
            [j*(n+1)+i for j in range(1,n) for i in range(1,n)]]
        reduced = LowRankTangent(tangent.local[active][:,active],
                                 tangent.vector[active],tangent.coefficient)
        solution, solve_audit = reduced.solve(force)
        full = np.zeros(reference.numel()); full[active]=solution
        observable = full[3*((n//2)*(n+1)+n//2)+2]; ref_value=PLATE_REFERENCE
    elif benchmark == "roof":
        template = roof_problem(n); active=template["active"]
        load = np.zeros(reference.numel()); load[2::3] = -90.0*(25*2*25*math.radians(40))*weights.reshape(-1)
        reduced=LowRankTangent(tangent.local[active][:,active],tangent.vector[active],tangent.coefficient)
        solution,solve_audit=reduced.solve(load[active]);full=np.zeros(reference.numel());full[active]=solution
        observable=-full[template["output_raw"]];ref_value=ROOF_REFERENCE
    else:
        template=m3_problem(n); template.update({"reference":reference,"hx":hx,"hy":hy,"targets":targets,"material":material})
        full,solve_audit=tangent.solve(template["force"],constraints=template["constraints"])
        observable=-float(template["point_b"][1]@full);ref_value=M3_REFERENCE
    row={"benchmark":benchmark,"n":n,"observable":float(observable),"reference":ref_value,
         "relative_error":abs(float(observable)-ref_value)/ref_value,
         "raw_scalar_states":int(reference.numel()),"local_tangent_nnz":int(tangent.local.nnz),
         "augmented_nnz":solve_audit["augmented_nnz"],**quality(reference),
         "minimum_logical_cell_area":float(np.min(signed)),
         "maximum_logical_cell_area":float(np.max(signed)),
         "reference_residual_inf":float(np.max(np.abs(residual))),
         "status":"SUCCESS"}
    return row, (reference.numpy()+full.reshape(reference.shape))


def run():
    raw=OUT/"raw_data"; figures=OUT/"figures"; config=OUT/"config"; logs=OUT/"logs"
    for path in (raw,figures,config,logs): path.mkdir(parents=True,exist_ok=True)
    (config/"experiment.json").write_text(json.dumps({"seed":SEED,"levels":LEVELS,"meshes":MESHES,
        "scope":"surface-preserving redistribution on identical plane/cylinder/sphere",
        "parameter_spacing":"nominal constant hx,hy retained; arbitrary variable-step derivation not claimed"},indent=2)+"\n")
    rows=[]
    for benchmark,n in MESHES.items():
        for kind in ("smooth","random"):
            for delta in LEVELS:
                xi,eta=logical_grid(n,delta,kind,benchmark)
                try:
                    row,deformed=solve_case(benchmark,n,xi,eta)
                    row.update({"method":"CB-Full","distortion_kind":kind,"delta":delta,"seed":SEED})
                    if delta==0.30:
                        reference,_,_=reference_for(benchmark,n,xi,eta)
                        np.savez_compressed(raw/f"{benchmark}_{kind}_030.npz",reference=reference.numpy(),deformed=deformed,xi=xi,eta=eta)
                except Exception as error:
                    row={"benchmark":benchmark,"n":n,"method":"CB-Full","distortion_kind":kind,
                         "delta":delta,"seed":SEED,"status":"FAILED","reason":f"{type(error).__name__}: {error}"}
                rows.append(row)
    keys=[]
    for row in rows: keys.extend(k for k in row if k not in keys)
    with (raw/"distortion.csv").open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=keys);w.writeheader();w.writerows(rows)
    (raw/"distortion.json").write_text(json.dumps(rows,indent=2,allow_nan=True)+"\n")
    markers={"plate":"o","roof":"s","M3":"^"}; styles={"smooth":"-","random":"--"}
    fig,ax=plt.subplots(figsize=(4.8,3.2))
    for b in MESHES:
        for kind in styles:
            sub=[r for r in rows if r["benchmark"]==b and r["distortion_kind"]==kind and r["status"]=="SUCCESS"]
            ax.plot([r["minimum_scaled_jacobian"] for r in sub],[r["relative_error"] for r in sub],
                    marker=markers[b],linestyle=styles[kind],label=f"{b}, {kind}")
    ax.set(xlabel="minimum scaled Jacobian",ylabel="relative displacement error",yscale="log")
    ax.grid(True,which="both",alpha=.25);ax.legend(fontsize=7,ncol=2)
    for ext in ("png","pdf","svg"): fig.savefig(figures/f"surface_preserving_distortion.{ext}",dpi=220,bbox_inches="tight")
    plt.close(fig)
    return rows


if __name__=="__main__":
    rows=run();print(json.dumps({"rows":len(rows),"failures":sum(r["status"]!="SUCCESS" for r in rows)},indent=2))

