"""Explicit analytical and synthetic morphology generators."""

from __future__ import annotations

from typing import Literal

import numpy as np
from scipy.ndimage import gaussian_filter

from .morphology import Morphology


def _coordinates(shape: tuple[int, int, int]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    axes = [np.arange(n, dtype=float) / n for n in shape]
    grid = np.meshgrid(*axes, indexing="ij")
    return grid[0], grid[1], grid[2]


def _quantile_phase(field: np.ndarray, fraction: float) -> np.ndarray:
    phase, _, _ = _rank_phase(field, fraction)
    return phase


def _rank_phase(field: np.ndarray, fraction: float) -> tuple[np.ndarray, float, int]:
    if not 0 < fraction < 1:
        raise ValueError("fraction must be strictly between zero and one")
    # Quantile ties can make the measured fraction drift. Deterministically select
    # the required number of tied voxels without changing the field definition.
    target = int(round(fraction * field.size))
    flat_order = np.argsort(field.ravel(), kind="stable")
    selected = flat_order[-target:]
    exact = np.zeros(field.size, dtype=bool)
    exact[selected] = True
    threshold = float(np.sort(field.ravel())[-target])
    return exact.reshape(field.shape), threshold, target


def lamellae(
    shape: tuple[int, int, int] = (48, 48, 48),
    fraction: float = 0.4,
    axis: int = 0,
    period_voxels: int | None = None,
) -> Morphology:
    """Generate planar conducting layers, periodic in all axes.

    The layer coordinate is a sawtooth and a deterministic rank threshold fixes
    the measured voxel fraction. A partial boundary plane can occur when the
    requested fraction is not exactly representable by whole voxel planes.
    """

    if axis not in (0, 1, 2):
        raise ValueError("axis must be 0, 1, or 2")
    n = shape[axis]
    period = period_voxels or n
    if period < 2 or period > n:
        raise ValueError("period_voxels must be between 2 and the axis length")
    coordinate_1d = np.arange(n, dtype=float) % period / period
    coordinate = np.broadcast_to(
        coordinate_1d.reshape(tuple(n if i == axis else 1 for i in range(3))), shape
    )
    phase, threshold, target = _rank_phase(-coordinate, fraction)
    return Morphology(
        phase=phase,
        name="lamellae",
        parameters={
            "fraction_target": fraction,
            "axis": axis,
            "period_voxels": period,
            "rank_field": "negative periodic layer coordinate",
            "rank_threshold": threshold,
            "target_voxels": target,
        },
    )


def cylinders(
    shape: tuple[int, int, int] = (48, 48, 48),
    fraction: float = 0.4,
    axis: Literal[0, 1, 2] = 2,
) -> Morphology:
    """Generate a square lattice of parallel, periodic cylinders.

    The cylinder axis is the transport phase's continuous direction. The
    transverse unit cell contains four identical cylinders, making the geometry
    explicit rather than claiming a generic block-copolymer cylinder phase.
    """

    if not 0 < fraction < 1:
        raise ValueError("fraction must be strictly between zero and one")
    coords = _coordinates(shape)
    transverse = [i for i in range(3) if i != axis]
    u, v = coords[transverse[0]], coords[transverse[1]]
    # Four cylinders per unit cell, centered at quarter and three-quarter points.
    centers = (0.25, 0.75)
    field = np.full(shape, np.inf)
    for cu in centers:
        for cv in centers:
            du = np.abs(u - cu)
            dv = np.abs(v - cv)
            du = np.minimum(du, 1.0 - du)
            dv = np.minimum(dv, 1.0 - dv)
            field = np.minimum(field, np.sqrt(du**2 + dv**2))
    phase, threshold, target = _rank_phase(-field, fraction)
    nominal_radius = np.sqrt(fraction / np.pi) / 2.0
    return Morphology(
        phase=phase,
        name="cylinders",
        parameters={
            "fraction_target": fraction,
            "axis": axis,
            "lattice": "square",
            "rank_field": "negative distance to periodic cylinder centers",
            "rank_threshold": threshold,
            "target_voxels": target,
            "nominal_continuum_radius_fraction_of_cell": nominal_radius,
        },
    )


def gyroid_level_set(
    shape: tuple[int, int, int] = (48, 48, 48),
    fraction: float = 0.4,
) -> Morphology:
    """Generate a thresholded triply-periodic gyroid level-set reference.

    This is a mathematical level-set fixture, not an SCFT or crystallization
    result. The standard first-harmonic field is thresholded to the requested
    measured volume fraction.
    """

    x, y, z = _coordinates(shape)
    two_pi = 2.0 * np.pi
    field = (
        np.sin(two_pi * x) * np.cos(two_pi * y)
        + np.sin(two_pi * y) * np.cos(two_pi * z)
        + np.sin(two_pi * z) * np.cos(two_pi * x)
    )
    phase, threshold, target = _rank_phase(field, fraction)
    return Morphology(
        phase=phase,
        name="gyroid_level_set",
        parameters={
            "fraction_target": fraction,
            "rank_field": "first-harmonic gyroid field",
            "rank_threshold": threshold,
            "target_voxels": target,
            "level_set_interpretation": "one thresholded level-set phase; not a double gyroid",
        },
    )


def bicontinuous_reference(
    shape: tuple[int, int, int] = (48, 48, 48),
    fraction: float = 0.4,
    seed: int = 20260910,
) -> Morphology:
    """Generate a periodic low-frequency Gaussian random-field fixture.

    The fixture is intended as a reproducible disordered reference. It is not
    labelled as a polymer morphology or a phase-field crystallization result.
    """

    rng = np.random.default_rng(seed)
    white = rng.standard_normal(shape)
    # Periodic convolution avoids a special boundary treatment at the edges.
    field = gaussian_filter(white, sigma=max(shape) / 14.0, mode="wrap")
    phase, threshold, target = _rank_phase(field, fraction)
    return Morphology(
        phase=phase,
        name="gaussian_random_field",
        parameters={
            "fraction_target": fraction,
            "seed": seed,
            "filter_mode": "periodic",
            "gaussian_sigma_voxels": max(shape) / 14.0,
            "rank_field": "periodic Gaussian-filtered white noise",
            "rank_threshold": threshold,
            "target_voxels": target,
        },
    )


def canonical_morphologies(
    shape: tuple[int, int, int] = (48, 48, 48), fraction: float = 0.4, seed: int = 20260910
) -> list[Morphology]:
    """Return the M1 comparison set."""

    return [
        lamellae(shape, fraction=fraction, axis=0),
        cylinders(shape, fraction=fraction, axis=2),
        gyroid_level_set(shape, fraction=fraction),
        bicontinuous_reference(shape, fraction=fraction, seed=seed),
    ]
