"""Periodic connected-component and percolation analysis."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .morphology import Morphology


class _UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = np.arange(size, dtype=np.int64)
        self.rank = np.zeros(size, dtype=np.int8)

    def find(self, item: int) -> int:
        root = item
        while self.parent[root] != root:
            root = int(self.parent[root])
        while self.parent[item] != item:
            next_item = int(self.parent[item])
            self.parent[item] = root
            item = next_item
        return root

    def union(self, first: int, second: int) -> None:
        left, right = self.find(first), self.find(second)
        if left == right:
            return
        if self.rank[left] < self.rank[right]:
            left, right = right, left
        self.parent[right] = left
        if self.rank[left] == self.rank[right]:
            self.rank[left] += 1


@dataclass(frozen=True)
class ConnectivityMetrics:
    components: int
    largest_component_fraction: float
    percolates: tuple[bool, bool, bool]
    phase_fraction: float

    def as_dict(self) -> dict[str, float | int | bool]:
        return {
            "components": self.components,
            "largest_component_fraction": self.largest_component_fraction,
            "percolates_x": self.percolates[0],
            "percolates_y": self.percolates[1],
            "percolates_z": self.percolates[2],
            "phase_fraction": self.phase_fraction,
        }


def analyze_connectivity(morphology: Morphology) -> ConnectivityMetrics:
    """Analyze face-connected transport phase with periodic edges where declared.

    A directional periodic percolation flag is true when one periodic connected
    component touches both opposing faces in that direction. For a nonperiodic
    axis, the same face-to-face test describes through-sample connectivity.
    """

    phase = morphology.phase
    shape = phase.shape
    active = np.flatnonzero(phase.ravel())
    if active.size == 0:
        return ConnectivityMetrics(0, 0.0, (False, False, False), 0.0)
    uf = _UnionFind(phase.size)
    for axis in range(3):
        current = np.argwhere(phase)
        neighbor = current.copy()
        valid = np.ones(len(current), dtype=bool)
        neighbor[:, axis] += 1
        if morphology.periodic_axes[axis]:
            neighbor[:, axis] %= shape[axis]
        else:
            valid = neighbor[:, axis] < shape[axis]
        current_ids = np.ravel_multi_index(current[valid].T, shape)
        neighbor_ids = np.ravel_multi_index(neighbor[valid].T, shape)
        valid_pairs = phase.ravel()[neighbor_ids]
        for left, right in zip(current_ids[valid_pairs], neighbor_ids[valid_pairs], strict=True):
            uf.union(int(left), int(right))

    roots = np.array([uf.find(int(item)) for item in active], dtype=np.int64)
    unique, counts = np.unique(roots, return_counts=True)
    component_count = int(unique.size)
    largest_fraction = float(counts.max() / active.size)
    percolates: list[bool] = []
    active_coords = np.argwhere(phase)
    active_roots = np.array([uf.find(int(np.ravel_multi_index(c, shape))) for c in active_coords])
    for axis in range(3):
        low = active_coords[:, axis] == 0
        high = active_coords[:, axis] == shape[axis] - 1
        percolates.append(bool(np.intersect1d(active_roots[low], active_roots[high]).size))
    return ConnectivityMetrics(
        components=component_count,
        largest_component_fraction=largest_fraction,
        percolates=tuple(percolates),
        phase_fraction=float(np.mean(phase)),
    )


def interfacial_area_density(morphology: Morphology) -> float:
    """Estimate binary interface area per physical volume from voxel faces."""

    phase = morphology.phase
    area = 0.0
    for axis, _spacing in enumerate(morphology.spacing):
        shifted = np.roll(phase, -1, axis=axis) if morphology.periodic_axes[axis] else None
        if shifted is None:
            slices_left = [slice(None)] * 3
            slices_right = [slice(None)] * 3
            slices_left[axis] = slice(0, -1)
            slices_right[axis] = slice(1, None)
            faces = phase[tuple(slices_left)] != phase[tuple(slices_right)]
        else:
            faces = phase != shifted
        face_area = np.prod([morphology.spacing[i] for i in range(3) if i != axis])
        area += float(np.count_nonzero(faces)) * face_area
    volume = float(np.prod(morphology.physical_size))
    return area / volume
