import numpy as np
import torch

from checkerboard_shell import (BoundaryConditions, CheckerboardMesh,
                                CheckerboardShell, IsotropicMaterial, NodalLoad)


def test_prescribed_sparse_tangent_and_reaction_balance():
    mesh = CheckerboardMesh.rectangle(width=1, height=.5, nx=4, ny=4)
    base = BoundaryConditions.clamped_boundary(mesh, "left")
    prescribed = torch.zeros(mesh.shape, dtype=torch.float64)
    prescribed[base.fixed_mask] = 1e-4
    bc = base.with_prescribed_displacement(prescribed)
    load = NodalLoad.total_force(mesh, (0, 0, -.01), mesh.boundaries["right"])
    model = CheckerboardShell(mesh, IsotropicMaterial(2e5, .3, .02), bc, load)
    q = bc.initial_coordinates(mesh.reference_vertices)
    tangent, sparse_r, _ = model.low_rank_tangent(q)
    dense_r = model.residual(q).detach().numpy()
    np.testing.assert_allclose(sparse_r, dense_r, atol=1e-8, rtol=1e-8)
    partial_load = 0.37
    _, sparse_partial, _ = model.low_rank_tangent(q, partial_load)
    dense_partial = model.residual(q, partial_load).detach().numpy()
    np.testing.assert_allclose(sparse_partial, dense_partial,
                               atol=1e-8, rtol=1e-8)
    np.testing.assert_allclose(sparse_partial-sparse_r,
                               (1-partial_load)*load.values[bc.free_mask].numpy(),
                               atol=1e-12, rtol=1e-12)
    direction = np.random.default_rng(0).normal(size=len(q))
    np.testing.assert_allclose(tangent.action(direction),
                               model.tangent(q).detach().numpy() @ direction,
                               atol=1e-5, rtol=1e-6)
    result = model.solve()
    assert result.status == "CONVERGED"
    residual, reaction = model.full_residual_and_reactions(result.coordinates)
    assert float(residual[bc.free_mask].abs().max()) < 1e-6
    np.testing.assert_allclose(reaction.sum((0, 1)).numpy(),
                               -load.values.sum((0, 1)).numpy(), atol=1e-6)


def test_patch_resultant_and_moment():
    mesh = CheckerboardMesh.rectangle(width=2, height=1, nx=5, ny=4)
    force = np.array([1.0, 2.0, 3.0])
    load = NodalLoad.quad_patch_force(mesh, 2, 1, .3, .7, force)
    xyz = mesh.reference_vertices.numpy().reshape(-1, 3)
    nodal = load.values.numpy().reshape(-1, 3)
    cell = np.array([xyz[1*(mesh.nx+1)+2], xyz[1*(mesh.nx+1)+3],
                     xyz[2*(mesh.nx+1)+3], xyz[2*(mesh.nx+1)+2]])
    weights = np.array([.7*.3, .3*.3, .3*.7, .7*.7])
    center = weights @ cell
    np.testing.assert_allclose(nodal.sum(0), force, atol=1e-14)
    np.testing.assert_allclose(np.cross(xyz, nodal).sum(0),
                               np.cross(center, force), atol=1e-14)


def test_dead_pressure_resultant_is_grid_independent():
    for nx,ny in ((4,4),(10,8)):
        mesh=CheckerboardMesh.rectangle(2,1,nx,ny)
        load=NodalLoad.dead_pressure(mesh,150.0)
        np.testing.assert_allclose(load.values.sum((0,1)).numpy(),
                                   [0,0,300],atol=1e-12)
        xyz=mesh.reference_vertices.numpy().reshape(-1,3)
        force=load.values.numpy().reshape(-1,3)
        np.testing.assert_allclose(np.cross(xyz,force).sum(0),
                                   [150,-300,0],atol=1e-12)


def test_prescribed_value_keeps_autograd_design_path():
    mesh=CheckerboardMesh.rectangle(1,.5,4,4)
    mask=BoundaryConditions.clamped_boundary(mesh,"left").fixed_mask.clone()
    mask[2,4,2]=True
    material=IsotropicMaterial(2e5,.3,.02)
    value=torch.tensor(1e-4,dtype=torch.float64,requires_grad=True)
    prescribed=torch.zeros(mesh.shape,dtype=torch.float64)
    prescribed[2,4,2]=value
    model=CheckerboardShell(mesh,material,BoundaryConditions(mask,prescribed))
    q=torch.zeros(int((~mask).sum()),dtype=torch.float64)
    grad,=torch.autograd.grad(model.potential(q),value)
    def energy(d):
        field=torch.zeros(mesh.shape,dtype=torch.float64)
        field[2,4,2]=d
        variant=CheckerboardShell(mesh,material,BoundaryConditions(mask,field))
        return float(variant.potential(q))
    step=1e-7
    numerical=(energy(1e-4+step)-energy(1e-4-step))/(2*step)
    np.testing.assert_allclose(float(grad),numerical,rtol=1e-6,atol=1e-8)
