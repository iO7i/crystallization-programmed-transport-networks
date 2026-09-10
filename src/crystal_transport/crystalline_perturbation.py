"""Controlled crystalline-like perturbations for Benchmark 002.

The functions in this module operate on an already-defined binary morphology. They
do not model crystal formation, kinetics, orientation, or chemistry. A selected
subset of the existing transport phase is only relabelled as a phenomenological
low-mobility subdomain.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.ndimage import distance_transform_edt

from .morphology import Morphology
from .transport import TransportResult

PlacementMode = Literal["random", "interface", "backbone", "low_criticality"]


@dataclass(frozen=True)
class PlacementResult:
    """Auditable selection of transport-phase voxels for a perturbation."""

    crystalline_mask: np.ndarray
    mode: PlacementMode
    seed: int
    requested_fraction_of_transport_phase: float
    eligible_voxel_count: int
    selected_voxel_count: int
    algorithm: str

    @property
    def actual_fraction_of_transport_phase(self) -> float:
        return self.selected_voxel_count / self.eligible_voxel_count

    @property
    def actual_total_material_fraction(self) -> float:
        return float(np.mean(self.crystalline_mask))

    @property
    def selection_digest(self) -> str:
        return hashlib.sha256(
            np.ascontiguousarray(self.crystalline_mask, dtype=np.uint8).tobytes()
        ).hexdigest()

    def as_dict(self) -> dict[str, object]:
        return {
            "placement_mode": self.mode,
            "seed": self.seed,
            "requested_crystalline_fraction_of_transport_phase": (
                self.requested_fraction_of_transport_phase
            ),
            "actual_crystalline_fraction_of_transport_phase": (
                self.actual_fraction_of_transport_phase
            ),
            "actual_crystalline_fraction_of_total_material": self.actual_total_material_fraction,
            "eligible_voxel_count": self.eligible_voxel_count,
            "selected_voxel_count": self.selected_voxel_count,
            "selection_algorithm": self.algorithm,
            "selection_digest": self.selection_digest,
        }


def build_diffusivity_field(
    morphology: Morphology,
    crystalline_mask: np.ndarray,
    d_mobile: float,
    d_crystal: float,
    d_matrix: float,
) -> np.ndarray:
    """Map matrix, mobile, and crystalline-like labels to a positive ``D(x)``."""

    mask = np.asarray(crystalline_mask, dtype=bool)
    if mask.shape != morphology.shape:
        raise ValueError("crystalline mask must have the morphology shape")
    if np.any(mask & ~morphology.phase):
        raise ValueError("crystalline-like voxels must be inside the transport phase")
    if min(d_mobile, d_crystal, d_matrix) <= 0:
        raise ValueError("all diffusivities must be positive")
    field = np.where(morphology.phase, d_mobile, d_matrix).astype(float)
    field[mask] = d_crystal
    return field


def _interface_distance(morphology: Morphology) -> np.ndarray:
    """Return periodic Euclidean distance to the nearest matrix voxel."""

    if np.all(morphology.phase):
        return np.full(morphology.shape, np.inf, dtype=float)
    tiled = np.tile(morphology.phase, (3, 3, 3))
    distances = distance_transform_edt(tiled, sampling=morphology.spacing)
    slices = tuple(slice(size, 2 * size) for size in morphology.shape)
    return distances[slices]


def baseline_flux_score(
    morphology: Morphology,
    baseline_results: tuple[TransportResult, TransportResult, TransportResult],
    d_mobile: float,
    d_matrix: float,
) -> np.ndarray:
    """Compute a non-circular baseline local flux-magnitude score.

    The score combines the three unperturbed scalar solves. It is used only to define
    deterministic high-flux (backbone) and low-flux placement hypotheses; the
    perturbed solve is never used to choose its own voxels.
    """

    local_d = np.where(morphology.phase, d_mobile, d_matrix).astype(float)
    squared = np.zeros(morphology.shape, dtype=float)
    for axis, result in enumerate(baseline_results):
        gradients: list[np.ndarray] = []
        for direction, spacing in enumerate(morphology.spacing):
            if direction == axis:
                gradient = np.gradient(
                    result.concentration,
                    spacing,
                    axis=direction,
                    edge_order=1,
                )
            else:
                gradient = (
                    np.roll(result.concentration, -1, axis=direction)
                    - np.roll(result.concentration, 1, axis=direction)
                ) / (2.0 * spacing)
            gradients.append(gradient)
        current = local_d * gradients[axis]
        squared += current**2
    return np.sqrt(squared)


def select_crystalline_like(
    morphology: Morphology,
    fraction_of_transport_phase: float,
    mode: PlacementMode,
    seed: int,
    baseline_flux: np.ndarray | None = None,
) -> PlacementResult:
    """Select an exact-count perturbation mask using a declared placement rule."""

    if not 0.0 <= fraction_of_transport_phase <= 1.0:
        raise ValueError("crystalline-like fraction must be between zero and one")
    if mode not in {"random", "interface", "backbone", "low_criticality"}:
        raise ValueError(f"unsupported placement mode: {mode}")
    eligible = np.flatnonzero(morphology.phase.ravel())
    selected_count = int(round(fraction_of_transport_phase * eligible.size))
    rng = np.random.default_rng(seed)
    if mode == "random":
        selected = rng.choice(eligible, size=selected_count, replace=False)
        algorithm = "seeded uniform sampling without replacement from transport phase"
    else:
        tie_break = rng.random(eligible.size)
        if mode == "interface":
            score = _interface_distance(morphology).ravel()[eligible]
            order = np.lexsort((tie_break, score))
            algorithm = "ascending periodic Euclidean distance to matrix interface"
        else:
            if baseline_flux is None or baseline_flux.shape != morphology.shape:
                raise ValueError("backbone and low-criticality modes require baseline_flux")
            score = baseline_flux.ravel()[eligible]
            descending = mode == "backbone"
            order = np.lexsort((tie_break, -score if descending else score))
            algorithm = (
                "descending baseline local flux-magnitude rank"
                if descending
                else "ascending baseline local flux-magnitude rank"
            )
        selected = eligible[order[:selected_count]]
    mask = np.zeros(morphology.phase.size, dtype=bool)
    mask[selected] = True
    return PlacementResult(
        crystalline_mask=mask.reshape(morphology.shape),
        mode=mode,
        seed=seed,
        requested_fraction_of_transport_phase=fraction_of_transport_phase,
        eligible_voxel_count=int(eligible.size),
        selected_voxel_count=selected_count,
        algorithm=algorithm,
    )
