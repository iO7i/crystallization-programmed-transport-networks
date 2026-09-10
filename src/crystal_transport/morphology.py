"""Canonical voxel morphology representation.

The project uses cell-centered, axis-ordered arrays with shape ``(nx, ny, nz)``.
Boolean ``phase`` values identify the phase assigned the transport coefficient
``D_transport`` in the benchmark.  This representation is intentionally small:
it carries enough metadata to make a generated morphology auditable without
introducing a general schema framework.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class Morphology:
    """A reproducible 3-D binary morphology on a Cartesian voxel grid."""

    phase: np.ndarray
    spacing: tuple[float, float, float] = (1.0, 1.0, 1.0)
    periodic_axes: tuple[bool, bool, bool] = (True, True, True)
    transport_phase: int = 1
    name: str = "unnamed"
    parameters: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        phase = np.asarray(self.phase, dtype=bool)
        if phase.ndim != 3:
            raise ValueError("phase must be a three-dimensional array")
        if any(size < 2 for size in phase.shape):
            raise ValueError("each morphology dimension must contain at least two voxels")
        if len(self.spacing) != 3 or any(value <= 0 for value in self.spacing):
            raise ValueError("spacing must contain three positive values")
        if len(self.periodic_axes) != 3:
            raise ValueError("periodic_axes must contain three booleans")
        object.__setattr__(self, "phase", phase)

    @property
    def shape(self) -> tuple[int, int, int]:
        return (int(self.phase.shape[0]), int(self.phase.shape[1]), int(self.phase.shape[2]))

    @property
    def volume_fraction(self) -> float:
        """Fraction of voxels assigned to the transport phase."""

        return float(np.mean(self.phase))

    @property
    def physical_size(self) -> tuple[float, float, float]:
        return (
            self.shape[0] * self.spacing[0],
            self.shape[1] * self.spacing[1],
            self.shape[2] * self.spacing[2],
        )

    def content_digest(self) -> str:
        """Return a stable digest of data and provenance-defining metadata."""

        metadata = {
            "spacing": self.spacing,
            "periodic_axes": self.periodic_axes,
            "transport_phase": self.transport_phase,
            "name": self.name,
            "parameters": self.parameters,
        }
        payload = json.dumps(metadata, sort_keys=True, separators=(",", ":")).encode()
        digest = hashlib.sha256()
        digest.update(payload)
        digest.update(np.ascontiguousarray(self.phase, dtype=np.uint8).tobytes())
        return digest.hexdigest()

    def metrics_metadata(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "shape": list(self.shape),
            "spacing": list(self.spacing),
            "periodic_axes": list(self.periodic_axes),
            "transport_phase": self.transport_phase,
            "volume_fraction": self.volume_fraction,
            "content_digest": self.content_digest(),
            "parameters": self.parameters,
        }

    def save(self, path: str | Path) -> None:
        """Write a compressed morphology file with JSON metadata."""

        target = Path(path)
        metadata = json.dumps(self.metrics_metadata(), sort_keys=True)
        np.savez_compressed(target, phase=self.phase.astype(np.uint8), metadata=metadata)

    @classmethod
    def load(cls, path: str | Path) -> Morphology:
        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata"]))
            return cls(
                phase=data["phase"].astype(bool),
                spacing=tuple(metadata["spacing"]),
                periodic_axes=tuple(metadata["periodic_axes"]),
                transport_phase=int(metadata["transport_phase"]),
                name=str(metadata["name"]),
                parameters=dict(metadata.get("parameters", {})),
            )
