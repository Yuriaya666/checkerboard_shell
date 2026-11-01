"""Convergence metrics and iteration records independent of a shell model."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any

import torch


def residual_norms(residual: torch.Tensor) -> tuple[float, float]:
    """Return infinity and Euclidean norms as Python floats."""
    value = residual.detach()
    if value.numel() == 0:
        return 0.0, 0.0
    return float(value.abs().max()), float(torch.linalg.vector_norm(value))


def convergence_metrics(residual: torch.Tensor, scale: float = 1.0) -> dict[str, float]:
    """Compute absolute and scale-normalized residual metrics."""
    residual_inf, residual_l2 = residual_norms(residual)
    scale = max(abs(float(scale)), 1.0e-30)
    return {
        "residual_inf": residual_inf,
        "residual_l2": residual_l2,
        "relative_residual_inf": residual_inf / scale,
        "relative_residual_l2": residual_l2 / scale,
    }


def is_converged(metrics: dict[str, float], absolute_tolerance: float,
                 relative_tolerance: float | None = None) -> bool:
    """Apply an absolute criterion, optionally combined with a relative one."""
    absolute = metrics["residual_inf"] <= float(absolute_tolerance)
    relative = (relative_tolerance is not None
                and metrics["relative_residual_inf"] <= float(relative_tolerance))
    if relative_tolerance is None:
        return absolute
    return absolute or relative


@dataclass
class ConvergenceMonitor:
    """Append-only monitor used by correctors and continuation drivers."""

    records: list[dict[str, Any]] = field(default_factory=list)

    def add(self, *, iteration: int, residual: torch.Tensor | None = None,
            metrics: dict[str, float] | None = None, **values: Any) -> dict[str, Any]:
        if metrics is None:
            if residual is None:
                raise ValueError("either residual or metrics is required")
            metrics = convergence_metrics(residual)
        record = {"iteration": int(iteration), **metrics, **values}
        self.records.append(record)
        return record

    @property
    def last(self) -> dict[str, Any] | None:
        return self.records[-1] if self.records else None

    @property
    def residual_inf(self) -> float:
        return math.inf if self.last is None else float(self.last["residual_inf"])
