"""Frozen position-only B/W Checkerboard geometry."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch


@dataclass
class TorchGeometry:
    t1_b: torch.Tensor
    t2_b: torch.Tensor
    normals_b: torch.Tensor
    first_form: torch.Tensor
    t1_w: torch.Tensor
    t2_w: torch.Tensor
    d1_normal: torch.Tensor
    d2_normal: torch.Tensor
    b_tilde: torch.Tensor
    second_form: torch.Tensor
    skew_diagnostic: torch.Tensor
    chi: torch.Tensor


@dataclass(frozen=True)
class EdgeMidpoints:
    """Physical horizontal/vertical edge midpoints of the raw grid."""

    horizontal: torch.Tensor
    vertical: torch.Tensor


def edge_midpoints(vertices: torch.Tensor) -> EdgeMidpoints:
    """Construct the physical midpoint representation used by B faces."""
    return EdgeMidpoints(
        horizontal=0.5 * (vertices[:, :-1] + vertices[:, 1:]),
        vertical=0.5 * (vertices[:-1] + vertices[1:]),
    )


def _normalize(values):
    return values / torch.linalg.vector_norm(values, dim=-1, keepdim=True)


def torch_geometry(vertices, hx, hy):
    """Evaluate the frozen B/W geometry on a structured vertex array."""
    q00 = vertices[:-1, :-1]
    q10 = vertices[:-1, 1:]
    q11 = vertices[1:, 1:]
    q01 = vertices[1:, :-1]
    ell = float(np.hypot(hx, hy))
    t1_b = (q11 - q00) / ell
    t2_b = (q01 - q10) / ell
    normals_b = _normalize(torch.linalg.cross(t1_b, t2_b, dim=-1))
    first_form = torch.stack((
        torch.stack(((t1_b * t1_b).sum(-1),
                     (t1_b * t2_b).sum(-1)), dim=-1),
        torch.stack(((t2_b * t1_b).sum(-1),
                     (t2_b * t2_b).sum(-1)), dim=-1),
    ), dim=-2)

    ne, nw = normals_b[1:, 1:], normals_b[1:, :-1]
    sw, se = normals_b[:-1, :-1], normals_b[:-1, 1:]
    t1_w = 0.5 * (t1_b[1:, 1:] + t1_b[:-1, :-1])
    t2_w = 0.5 * (t2_b[1:, :-1] + t2_b[:-1, 1:])
    d1_normal = (ne - sw) / ell
    d2_normal = (nw - se) / ell
    b_tilde = -torch.stack((
        torch.stack(((d1_normal * t1_w).sum(-1),
                     (d1_normal * t2_w).sum(-1)), dim=-1),
        torch.stack(((d2_normal * t1_w).sum(-1),
                     (d2_normal * t2_w).sum(-1)), dim=-1),
    ), dim=-2)
    second_form = 0.5 * (b_tilde + b_tilde.transpose(-1, -2))
    skew = 0.5 * (b_tilde - b_tilde.transpose(-1, -2))
    chi = ne - nw + sw - se
    return TorchGeometry(
        t1_b=t1_b, t2_b=t2_b, normals_b=normals_b,
        first_form=first_form, t1_w=t1_w, t2_w=t2_w,
        d1_normal=d1_normal, d2_normal=d2_normal,
        b_tilde=b_tilde, second_form=second_form,
        skew_diagnostic=skew, chi=chi)

