"""Comparison contract for independently validated shell backends.

No MidedgeTan/MITC4 surrogate is supplied.  An adapter must prove a matching
mesh, load, observable and constrained physical subspace before reporting data.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np

from common import BenchmarkCase


@dataclass(frozen=True)
class BaselineResult:
    method: str
    status: str
    observable: float
    total_scalar_states: int
    active_scalar_states: int
    tangent_nnz: int | None
    tangent_density: float | None
    assembly_seconds: float | None
    solve_seconds: float | None
    peak_memory_bytes: int | None
    linear_iterations: int | None
    provenance: str


class VerifiedShellAdapter(Protocol):
    """Only a separately tested implementation should implement this protocol."""

    name: str
    provenance: str

    def solve_reference_linearized(self, case: BenchmarkCase) -> BaselineResult:
        """Use the case's exact geometry, physical BC, dead loads and observable."""

    def spectral_reference_tangent(
        self, case: BenchmarkCase, count: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return the smallest eigenpairs after the same rigid/gauge reduction."""


AVAILABLE_BASELINES: dict[str, VerifiedShellAdapter] = {}


def require_verified_baseline(name: str) -> VerifiedShellAdapter:
    if name not in AVAILABLE_BASELINES:
        raise RuntimeError(
            f"{name} has no verified in-project or dependency-provided adapter; "
            "fabricated comparison values are forbidden")
    return AVAILABLE_BASELINES[name]

