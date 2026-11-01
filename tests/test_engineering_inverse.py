import numpy as np
import torch

from checkerboard_shell import (BoundaryConditions, CheckerboardMesh,
                                CheckerboardShell, IsotropicMaterial)
from experiments.engineering_inverse.common import positional_beam_stiffness
from experiments.engineering_inverse.automotive_panel import build_case
from experiments.engineering_inverse.mount_springs import (
    patch_weights, bracket_matrices, checkerboard_forward)
from experiments.engineering_inverse.automotive_spring_design import response_sensitivity
from experiments.engineering_inverse.aircraft_skin import (
    build_case as build_aircraft, layout_operator, prepare_layout_compliance)
from experiments.engineering_inverse.common import checkerboard_displacement_operator


def test_position_beam_rigid_motion_and_cantilever_scaling():
    xyz=np.zeros((9,3)); xyz[:,0]=np.linspace(0,1,9)
    E,A,I=70e9,1e-5,1e-10
    K=positional_beam_stiffness(xyz,range(9),E,A,I)
    translation=np.tile([0,0,1.0],9)
    np.testing.assert_allclose(K@translation,0,atol=1e-8)
    free=np.array([3*n+2 for n in range(2,9)])
    rhs=np.zeros(len(free));rhs[-1]=1.0
    tip=np.linalg.solve(K[np.ix_(free,free)],rhs)[-1]
    assert tip>0
    # The first two fixed nodes approximate a clamped slope. The beam
    # response must be close to the Euler-Bernoulli cantilever value.
    assert abs(tip-(1.0/(3*E*I)))/(1.0/(3*E*I))<.35


def test_implicit_displacement_derivative_matches_central_difference():
    mesh=CheckerboardMesh.rectangle(1,.5,4,4)
    material=IsotropicMaterial(2e5,.3,.02)
    mask=BoundaryConditions.clamped_boundary(mesh,"left").fixed_mask.clone()
    controlled=(2,3,2)
    mask[controlled]=True
    all_free=CheckerboardShell(mesh,material,
        BoundaryConditions(torch.zeros(mesh.shape,dtype=torch.bool)))
    K,_,_=all_free.low_rank_tangent(torch.zeros(3*mesh.num_vertices,dtype=torch.float64))
    full=np.asarray(K.local+K.coefficient*np.outer(K.vector,K.vector))
    fixed=np.flatnonzero(mask.numpy().reshape(-1))
    free=np.flatnonzero((~mask).numpy().reshape(-1))
    dof=np.ravel_multi_index(controlled,mesh.shape)
    response_dof=np.ravel_multi_index((3,4,2),mesh.shape)
    row=int(np.where(free==response_dof)[0][0])
    analytic=-np.linalg.solve(full[np.ix_(free,free)],full[free,dof])[row]
    def response(d):
        prescribed=torch.zeros(mesh.shape,dtype=torch.float64)
        prescribed[controlled]=d
        model=CheckerboardShell(mesh,material,BoundaryConditions(mask,prescribed))
        result=model.solve()
        assert result.status=="CONVERGED"
        return float((result.deformed_vertices-mesh.reference_vertices).reshape(-1)[response_dof])
    step=1e-6
    numerical=(response(step)-response(-step))/(2*step)
    np.testing.assert_allclose(analytic,numerical,atol=2e-3,rtol=2e-3)


def test_mount_regions_are_physical_and_distinct():
    for nx,ny in ((12,8),(20,16)):
        case=build_case(nx,ny)
        patches=case["mount_patches"]
        assert len(patches)==4
        assert all(patches)
        assert len(set().union(*map(set,patches)))==sum(map(len,patches))
        xyz=case["mesh"].reference_vertices.numpy().reshape(-1,3)
        centers=np.array([xyz[patch,:2].mean(axis=0) for patch in patches])
        np.testing.assert_allclose(centers[:,0],
            np.array([.15,.40,.65,.85])*.8,atol=.5*.8/nx)
        material_y=np.array([np.mean(np.asarray(patch)//(nx+1))*.45/ny
                             for patch in patches])
        np.testing.assert_allclose(material_y,.75*.45,atol=.5*.45/ny)


def test_area_integrated_bracket_weights_preserve_area_and_centroid():
    for nx,ny in ((12,8),(20,16)):
        case=build_case(nx,ny)
        W=patch_weights(case)
        np.testing.assert_allclose(W.sum(axis=1),1,atol=1e-12)
        x=np.tile(np.linspace(0,case["width_m"],nx+1),ny+1)
        y=np.repeat(np.linspace(0,case["height_m"],ny+1),nx+1)
        np.testing.assert_allclose(W@x,np.array([.15,.40,.65,.85])*.8,
                                   atol=1e-12)
        np.testing.assert_allclose(W@y,np.full(4,.75*.45),atol=1e-12)


def test_aircraft_shared_compliance_matches_direct_elimination():
    case=build_aircraft(12,24)
    base=checkerboard_displacement_operator(case["mesh"],case["material"],
        case["spar"],[],case["O"],case["load"],beam_lines=case["beam_lines"])
    fast=layout_operator(case,(0,2,4),prepare_layout_compliance(case,base))
    direct=layout_operator(case,(0,2,4),
                           {"K":base["K"],"residual":base["residual"]})
    for name in ("q0","Q","response_map","actuator_reaction0",
                 "actuator_reaction_map"):
        np.testing.assert_allclose(fast[name],direct[name],rtol=1e-7,atol=1e-7)


def test_aircraft_sparse_rank_one_shared_matches_dense_reference():
    from experiments.engineering_inverse.aircraft_sparse import sparse_base_operator
    case=build_aircraft(12,24)
    dense=checkerboard_displacement_operator(case["mesh"],case["material"],
        case["spar"],[],case["O"],case["load"],beam_lines=case["beam_lines"])
    sparse=sparse_base_operator(case)
    a=layout_operator(case,(0,2,4),prepare_layout_compliance(case,dense))
    b=layout_operator(case,(0,2,4),prepare_layout_compliance(case,sparse))
    for name in ("response0","response_map","actuator_reaction0",
                 "actuator_reaction_map"):
        np.testing.assert_allclose(a[name],b[name],rtol=1e-8,atol=1e-8)


def test_mitc4_shared_factor_matches_direct_aircraft_layout():
    from experiments.engineering_inverse.mitc4_shared import build_shared,layout_operator as layout_mt
    from experiments.engineering_inverse.efficiency_aircraft import _mitc4_direct_layout
    case=build_aircraft(12,24)
    shared=build_shared(case)
    for subset in ((0,),(0,2,4)):
        direct=_mitc4_direct_layout(case,subset)
        fast=layout_mt(case,shared,subset)
        for name in ("response0","response_map","actuator_reaction0",
                     "actuator_reaction_map"):
            np.testing.assert_allclose(direct[name],fast[name],rtol=1e-7,atol=1e-7)


def test_bracket_implicit_shape_and_reaction_derivatives():
    case=build_case(12,8)
    controls=np.array([-.09,-.295,-.106,-.257,-.161,.366,-.121,.216])*1e-3
    baseline=checkerboard_forward(case,controls)
    assert baseline["status"]=="CONVERGED"
    K,D,W=bracket_matrices(case)
    J,dF=response_sensitivity(case,baseline,K,D,W)
    step=1e-5
    perturb=np.zeros(8);perturb[0]=step
    initial=torch.as_tensor(baseline["displacement"]).reshape(case["mesh"].shape)
    plus=checkerboard_forward(case,controls+perturb,initial=initial)
    minus=checkerboard_forward(case,controls-perturb,initial=initial)
    assert plus["status"]==minus["status"]=="CONVERGED"
    edge_fd=case["O"]@(plus["displacement"]-minus["displacement"]).ravel()/(2*step)
    force_fd=(plus["mount_reaction"]-minus["mount_reaction"])/(2*step)
    np.testing.assert_allclose(J[:,0],edge_fd,rtol=1e-3,atol=3e-4)
    np.testing.assert_allclose(dF[:,:,0],force_fd,rtol=1e-3,atol=.1)


def test_automotive_opposite_sublattice_translation_is_gauged():
    case=build_case(12,8)
    assert case["gauge_x_dof"]//3 not in case["datum"]
    initial=torch.zeros(case["mesh"].shape,dtype=torch.float64)
    for j in range(9):
        for i in range(13):
            if (i+j)%2:
                initial[j,i,0]=.001
    zero=checkerboard_forward(case,np.zeros(8))
    shifted=checkerboard_forward(case,np.zeros(8),initial=initial)
    assert zero["status"]==shifted["status"]=="CONVERGED"
    np.testing.assert_allclose(zero["displacement"],shifted["displacement"],
                               rtol=0,atol=1e-8)
    assert abs(shifted["gauge_reaction_N"])<1e-5
