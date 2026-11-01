"""Independent corotational shell/frame check of the selected aircraft design.

ASDShellQ4 and corotational elasticBeamColumn use OpenSees geometry updates.
The beam/shell rotation coupling differs from the position-line model and is
reported explicitly; this is a cross-model reference, not an identical element.
"""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np

from experiments.cmame_round4.mitc4 import _ops
from experiments.engineering_inverse.aircraft_skin import OUT, build_case


def forward(case, selected, strokes, steps=8):
    ops=_ops();ops.wipe()
    ops.model("basic","-ndm",3,"-ndf",6)
    xyz=case["mesh"].reference_vertices.detach().numpy().reshape(-1,3)
    for node,point in enumerate(xyz,1):ops.node(node,*map(float,point))
    mat=case["material"]
    ops.section("ElasticMembranePlateSection",1,mat.young,mat.poisson,
                mat.thickness,0.0)
    for tag,quad in enumerate(case["mesh"].quads,1):
        ops.element("ASDShellQ4",tag,*[int(v)+1 for v in quad],1,
                    "-corotational")
    tag=len(case["mesh"].quads)+1
    for line_no,(line,E,A,I) in enumerate(case["beam_lines"],1):
        ops.geomTransf("Corotational",line_no,0.0,0.0,1.0)
        G=E/(2*(1+.3))
        for a,b in zip(line[:-1],line[1:]):
            ops.element("elasticBeamColumn",tag,int(a)+1,int(b)+1,
                        A,E,G,2*I,I,I,line_no)
            tag+=1
    spar=set(case["spar"])
    groups=[[3*n+2 for n in case["patches"][station]] for station in selected]
    design={d:float(value) for group,value in zip(groups,strokes) for d in group}
    for node in range(len(xyz)):
        fix=[int(3*node+c in spar and 3*node+c not in design)
             for c in range(3)]
        if any(fix):ops.fix(node+1,*fix,0,0,0)
    ops.timeSeries("Linear",1);ops.pattern("Plain",1,1)
    for dof,value in design.items():ops.sp(dof//3+1,dof%3+1,value)
    load=case["load"].values.detach().numpy().reshape(-1,3)
    for node,force in enumerate(load,1):
        if np.any(force):ops.load(node,*map(float,force),0,0,0)
    ops.system("SparseGeneral");ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.integrator("LoadControl",1/steps)
    ops.test("NormUnbalance",1e-6,80)
    ops.algorithm("NewtonLineSearch")
    ops.analysis("Static")
    status=ops.analyze(steps)
    if status!=0:
        raise RuntimeError(f"OpenSees corotational analysis failed: {status}")
    q=np.array([ops.nodeDisp(node+1,c+1) for node in range(len(xyz))
                for c in range(3)])
    ops.reactions()
    reaction=[sum(ops.nodeReaction(d//3+1,d%3+1) for d in group)
              for group in groups]
    return dict(status=int(status),displacement_m=q.tolist(),
                attachment_reaction_N=reaction,
                load_factor=float(ops.getTime()))


def run(meshes=((12,24),(20,40)),selected=(0,2,4)):
    source=json.loads((OUT/"aircraft_nonlinear_final_12x24.json").read_text())
    report=dict(reference="OpenSees ASDShellQ4 -corotational plus corotational frame beams",
                beam_connection_difference="shared rotational DOFs versus position-line slopes",
                stations=list(selected),meshes={})
    for nx,ny in meshes:
        case=build_case(nx,ny)
        rows=[]
        for target_id,old in enumerate(source["cases"]):
            strokes=np.array(old["strokes_m"])
            result=forward(case,selected,strokes)
            q=np.array(result["displacement_m"])
            error=case["O"]@q-case["targets"][target_id]
            result.update(rms_error_mm=float(np.sqrt(np.mean(error**2))*1e3),
                          max_error_mm=float(np.max(np.abs(error))*1e3),
                          max_displacement_mm=float(np.max(np.linalg.norm(q.reshape(-1,3),axis=1))*1e3))
            rows.append(result)
            print(nx,ny,target_id,result["status"],result["rms_error_mm"],
                  result["max_error_mm"],flush=True)
        report["meshes"][f"{nx}x{ny}"]=rows
        (OUT/"aircraft_corotational_reference.json").write_text(json.dumps(report,indent=2)+"\n")
    return report


if __name__=="__main__":run()
