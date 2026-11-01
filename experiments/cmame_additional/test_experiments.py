"""Focused regression tests for the independent experimental harness."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import torch

from common import (alternating_vector, benchmark_case, mesh_quality,
                    parameter_grid)
from checkerboard_shell.geometry.core import edge_midpoints


def test_reference_quad_quality_is_positive():
    for benchmark in ("plate", "roof", "m3"):
        quality = mesh_quality(benchmark_case(benchmark, 8).mesh)
        assert quality["inverted_corner_count"] == 0
        assert quality["minimum_scaled_jacobian"] > 0.9


def test_random_and_smooth_distortion_match_interior_rms():
    for delta in (0.05, 0.2, 0.4):
        base = parameter_grid(8, 8)
        smooth = parameter_grid(8, 8, distortion=delta, distortion_kind="smooth")
        random = parameter_grid(8, 8, distortion=delta, distortion_kind="random")
        for index in range(2):
            a = smooth[index][1:-1, 1:-1]-base[index][1:-1, 1:-1]
            b = random[index][1:-1, 1:-1]-base[index][1:-1, 1:-1]
            assert np.isclose(np.linalg.norm(a), np.linalg.norm(b), rtol=1e-13)


def test_load_and_observable_annihilate_gauge():
    for benchmark in ("plate", "roof", "m3"):
        case = benchmark_case(benchmark, 8, distortion=0.2,
                              distortion_kind="smooth")
        parity = alternating_vector(case.mesh.shape)
        for component in range(3):
            direction = torch.zeros(3, dtype=torch.float64)
            direction[component] = 1.0
            mode = parity*direction
            assert abs(float((case.load*mode).sum())) < 1e-10
            assert abs(float((case.observable_weights*mode).sum())) < 1e-12


def test_edge_midpoints_are_exactly_raw_gauge_invariant():
    case = benchmark_case("m3", 8)
    q = case.mesh.reference_vertices
    direction = torch.tensor([1.0, 0.7, -0.4], dtype=q.dtype)
    shifted = q+0.25*alternating_vector(q.shape)*direction
    a, b = edge_midpoints(q), edge_midpoints(shifted)
    assert torch.allclose(a.horizontal, b.horizontal, atol=1e-14)
    assert torch.allclose(a.vertical, b.vertical, atol=1e-14)

