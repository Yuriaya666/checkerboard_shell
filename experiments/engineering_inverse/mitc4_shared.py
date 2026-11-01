"""OpenSees ShellMITC4 shared-factor aircraft layout operator.

One spar-clamped OpenSees model contains a pressure load case followed by
individual candidate-node unit-load cases in Path time series.  The domain
topology does not change between steps, so Linear -factorOnce reuses its
factorization.  Candidate compliance then enforces any layout through a
small dense Lagrange-multiplier solve, exactly as the Checkerboard path does.
"""
from __future__ import annotations

import time
import numpy as np

from experiments.cmame_round4.mitc4 import _ops


def _pattern(ops, tag, times, active_index, node, vector):
    values=[float(k==active_index) for k in range(len(times))]
    ops.timeSeries("Path",tag,"-time",*times,"-values",*values)
    ops.pattern("Plain",tag,tag)
    ops.load(int(node),*map(float,vector))


def build_shared(case):
    """Return baseline and candidate compliance using one sparse factorization."""
    ops=_ops();ops.wipe()
    start=time.perf_counter()
    ops.model("basic","-ndm",3,"-ndf",6)
    xyz=case["mesh"].reference_vertices.detach().numpy().reshape(-1,3)
    n=len(xyz)
    for node,point in enumerate(xyz,1):ops.node(node,*map(float,point))
    material=case["material"]
    ops.section("ElasticMembranePlateSection",1,material.young,
                material.poisson,material.thickness,0.0)
    for tag,quad in enumerate(case["mesh"].quads,1):
        ops.element("ShellMITC4",tag,*[int(v)+1 for v in quad],1)
    next_tag=len(case["mesh"].quads)+1
    for line_number,(line,young,area,inertia) in enumerate(case["beam_lines"],1):
        ops.geomTransf("Linear",line_number,0.0,0.0,1.0)
        shear=young/(2*(1+.3))
        for a,b in zip(line[:-1],line[1:]):
            ops.element("elasticBeamColumn",next_tag,int(a)+1,int(b)+1,
                        area,young,shear,2*inertia,inertia,inertia,line_number)
            next_tag+=1
    spar=set(case["spar"])
    for node in range(n):
        fix=[int(3*node+c in spar) for c in range(3)]
        if any(fix):ops.fix(node+1,*fix,0,0,0)
    candidate=np.array(sorted({int(node) for patch in case["patches"] for node in patch}))
    steps=len(candidate)+1
    times=list(map(float,range(1,steps+1)))
    load=case["load"].values.detach().numpy().reshape(-1,3)
    tag=1
    for node,force in enumerate(load,1):
        if np.any(force):
            _pattern(ops,tag,times,0,node,(*force,0,0,0));tag+=1
    for column,node in enumerate(candidate):
        _pattern(ops,tag,times,column+1,node+1,(0,0,1,0,0,0));tag+=1
    construction_seconds=time.perf_counter()-start
    ops.system("SparseGeneral");ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.integrator("LoadControl",1.0)
    ops.algorithm("Linear","-factorOnce")
    ops.analysis("Static")
    observations=case["O"]
    observed_nodes=np.flatnonzero(observations.sum(axis=0).reshape(-1,3)[:,2])
    responses=np.zeros((observations.shape[0],steps))
    candidate_displacement=np.zeros((len(candidate),steps))
    counts=[];step_seconds=[]
    for step in range(steps):
        started=time.perf_counter()
        status=ops.analyze(1)
        step_seconds.append(time.perf_counter()-started)
        if status!=0:raise RuntimeError(f"OpenSees shared factor step {step} failed: {status}")
        counts.append(int(ops.numFact()))
        responses[:,step]=[ops.nodeDisp(int(node)+1,3) for node in observed_nodes]
        candidate_displacement[:,step]=[
            ops.nodeDisp(int(node)+1,3) for node in candidate]
    if len(set(counts))!=1:
        raise RuntimeError(f"OpenSees refactorized despite -factorOnce: {counts}")
    # SparseGeneral reports zero from numFact() in this OpenSeesPy build.
    # The -factorOnce contract is documented, but that counter cannot audit it.
    return dict(candidate=candidate,candidate_lookup={int(node):k for k,node in enumerate(candidate)},
                response_baseline=responses[:,0],response_influence=responses[:,1:],
                candidate_baseline=candidate_displacement[:,0],
                compliance=candidate_displacement[:,1:],
                construction_seconds=construction_seconds,
                first_analysis_seconds=step_seconds[0],
                candidate_response_seconds=sum(step_seconds[1:]),
                step_seconds=step_seconds,factorization_count_reported=counts[-1],
                total_dofs=6*n,free_dofs=6*n-len(spar),
                candidate_count=len(candidate))


def layout_operator(case,shared,selected):
    """Evaluate one layout from the shared candidate compliance matrix."""
    groups=[[int(node) for node in case["patches"][station]] for station in selected]
    chosen=sorted(set(node for group in groups for node in group))
    indices=np.array([shared["candidate_lookup"][node] for node in chosen])
    C=shared["compliance"][np.ix_(indices,indices)]
    B=np.zeros((len(chosen),len(groups)))
    lookup={node:k for k,node in enumerate(chosen)}
    for col,group in enumerate(groups):B[[lookup[node] for node in group],col]=1
    lam0=np.linalg.solve(C,-shared["candidate_baseline"][indices])
    lammap=np.linalg.solve(C,B)
    H=shared["response_influence"][:,indices]
    reaction0=np.array([lam0[[lookup[node] for node in group]].sum()
                        for group in groups])
    reaction_map=np.array([lammap[[lookup[node] for node in group]].sum(axis=0)
                           for group in groups])
    return dict(response0=shared["response_baseline"]+H@lam0,
                response_map=H@lammap,
                actuator_reaction0=reaction0,
                actuator_reaction_map=reaction_map)
