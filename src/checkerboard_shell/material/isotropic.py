"""Isotropic plane-stress material coefficients used by the paper model."""

from __future__ import annotations

from dataclasses import dataclass

import torch


@dataclass(frozen=True)
class IsotropicMaterial:
    """Young modulus, Poisson ratio and shell thickness (SI units)."""

    young: float = 1.2e6
    poisson: float = 0.0
    thickness: float = 0.1

    def __post_init__(self):
        if self.young <= 0.0 or self.thickness <= 0.0:
            raise ValueError("young and thickness must be positive")
        if not (-1.0 < self.poisson < 1.0):
            raise ValueError("poisson must lie in (-1, 1)")

    @property
    def lame_alpha(self) -> float:
        return float(self.young * self.poisson / (1.0 - self.poisson**2))

    @property
    def lame_mu(self) -> float:
        return float(self.young / (2.0 * (1.0 + self.poisson)))


def lame_surface_coefficients(young: float, poisson: float) -> tuple[float, float]:
    if not (-1.0 < poisson < 1.0):
        raise ValueError("Poisson ratio must lie in (-1,1)")
    return (float(young * poisson / (1.0 - poisson * poisson)),
            float(young / (2.0 * (1.0 + poisson))))


def stvk_density(change, reference_metric, alpha, beta):
    """LibShell-compatible ``.5 alpha tr(M)^2 + beta tr(M^2)`` density."""
    inverse = torch.linalg.inv(reference_metric)
    mixed = inverse @ change
    trace = torch.diagonal(mixed, dim1=-2, dim2=-1).sum(-1)
    squared = torch.diagonal(mixed @ mixed, dim1=-2, dim2=-1).sum(-1)
    return 0.5*float(alpha)*trace.square() + float(beta)*squared
