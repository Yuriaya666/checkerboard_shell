"""Isolated translational-line versus six-DOF frame cantilever comparison."""
from __future__ import annotations

import json
import numpy as np

from experiments.engineering_inverse.aircraft_skin import OUT
from experiments.engineering_inverse.common import positional_beam_stiffness
from experiments.cmame_round4.mitc4 import _ops


def run(cell_counts=(12,16,20,24,40)):
    E=70e9;width=.008;thickness=.001;length=.18
    area=width*thickness;inertia=width*thickness**3/12
    rows=[]
    for cells in cell_counts:
        xyz=np.column_stack((np.linspace(0,length,cells+1),
                             np.zeros(cells+1),np.zeros(cells+1)))
        K=positional_beam_stiffness(xyz,range(cells+1),E,area,inertia)
        free=np.array([3*node+2 for node in range(2,cells+1)])
        load=np.zeros(len(free));load[-1]=1
        tip_position=float(np.linalg.solve(K[np.ix_(free,free)],load)[-1])
        ops=_ops();ops.wipe();ops.model("basic","-ndm",3,"-ndf",6)
        for node,point in enumerate(xyz,1):ops.node(node,*point)
        ops.geomTransf("Linear",1,0,0,1)
        shear=E/(2*(1+.3))
        for node in range(1,cells+1):
            ops.element("elasticBeamColumn",node,node,node+1,
                        area,E,shear,2*inertia,inertia,inertia,1)
        ops.fix(1,1,1,1,1,1,1)
        ops.timeSeries("Linear",1);ops.pattern("Plain",1,1)
        ops.load(cells+1,0,0,1,0,0,0)
        ops.system("SparseGeneral");ops.numberer("RCM")
        ops.constraints("Transformation");ops.integrator("LoadControl",1)
        ops.algorithm("Linear");ops.analysis("Static")
        if ops.analyze(1)!=0:raise RuntimeError("frame cantilever failed")
        tip_frame=float(ops.nodeDisp(cells+1,3))
        rows.append(dict(cells=cells,position_tip_m_per_N=tip_position,
                         frame_tip_m_per_N=tip_frame,
                         position_over_frame=tip_position/tip_frame))
    payload=dict(test="isolated 0.18 m chordwise cantilever, 1 N end load",
        section=dict(E_Pa=E,area_m2=area,I_m4=inertia),
        boundary_difference="position beam clamps first two translations; frame clamps first node translations and rotations",
        rows=rows)
    (OUT/"beam_kinematics_audit.json").write_text(json.dumps(payload,indent=2)+"\n")
    return payload


if __name__=="__main__":print(run()["rows"])
