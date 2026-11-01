import numpy as np
import torch

from checkerboard_shell.material.isotropic import IsotropicMaterial
from checkerboard_shell.mechanics.energy import make_targets
from checkerboard_shell.mechanics.low_rank import assemble_paper_low_rank_tangent


def test_four_completion_variants_add_exactly_at_deformed_state():
    n = 4
    x = torch.linspace(0, 1, n+1, dtype=torch.float64)
    yy, xx = torch.meshgrid(x, x, indexing="ij")
    reference = torch.stack((xx, yy, 0.1*xx*yy), -1)
    displacement = torch.stack((0.01*xx*yy, -0.005*xx**2,
                                0.02*torch.sin(xx+yy)), -1)
    targets = make_targets(reference, 1/n, 1/n)
    material = IsotropicMaterial(3.0, 0.25, 0.02)
    tangents, gradients = {}, {}
    for variant in ("CB-0", "CB-M", "CB-W", "CB-Full"):
        tangents[variant], gradients[variant], _ = assemble_paper_low_rank_tangent(
            displacement, reference, 1/n, 1/n, targets, material, variant=variant)
    direction = np.random.default_rng(17).normal(size=reference.numel())
    expected = (tangents["CB-M"].action(direction)
                + tangents["CB-W"].action(direction)
                - tangents["CB-0"].action(direction))
    np.testing.assert_allclose(tangents["CB-Full"].action(direction), expected,
                               rtol=2e-7, atol=2e-8)
    np.testing.assert_allclose(gradients["CB-Full"],
                               gradients["CB-M"] + gradients["CB-W"] - gradients["CB-0"],
                               rtol=2e-7, atol=2e-8)


def test_original_curved_benchmark_responses_are_unchanged():
    from experiments.cmame_round2.recovered_drivers import solve_m3, solve_roof

    np.testing.assert_allclose(solve_roof(12)["displacement"],
                               0.33703867932120235, rtol=2e-9)
    np.testing.assert_allclose(solve_m3(8)["displacement_raw"],
                               0.07779532294558913, rtol=2e-9)
