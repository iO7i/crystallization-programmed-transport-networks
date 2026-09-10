"""Cell-centred finite-volume steady diffusion solver."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import cg

from .morphology import Morphology


@dataclass(frozen=True)
class TransportResult:
    axis: int
    effective_diffusivity: float
    flux: float
    residual_norm: float
    iterations: int
    concentration: np.ndarray


def _harmonic_mean(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    denominator = left + right
    return np.divide(
        2.0 * left * right,
        denominator,
        out=np.zeros_like(denominator),
        where=denominator > 0,
    )


def _edge_ids(
    shape: tuple[int, int, int], axis: int, periodic: bool
) -> tuple[np.ndarray, np.ndarray]:
    ids = np.arange(np.prod(shape), dtype=np.int64).reshape(shape)
    if periodic:
        return ids.ravel(), np.roll(ids, -1, axis=axis).ravel()
    left = [slice(None)] * 3
    right = [slice(None)] * 3
    left[axis] = slice(0, -1)
    right[axis] = slice(1, None)
    return ids[tuple(left)].ravel(), ids[tuple(right)].ravel()


def _assemble_system(
    morphology: Morphology,
    d_transport: float,
    d_matrix: float,
    axis: int,
    low_value: float,
    high_value: float,
) -> tuple[csr_matrix, np.ndarray, np.ndarray]:
    shape = morphology.shape
    spacing = np.asarray(morphology.spacing, dtype=float)
    diffusivity = np.where(morphology.phase, d_transport, d_matrix).astype(float)
    diffusivity_flat = diffusivity.ravel()
    row_parts: list[np.ndarray] = []
    col_parts: list[np.ndarray] = []
    data_parts: list[np.ndarray] = []
    diagonal = np.zeros(diffusivity.size, dtype=float)
    rhs = np.zeros(diffusivity.size, dtype=float)

    for direction in range(3):
        left_ids, right_ids = _edge_ids(shape, direction, morphology.periodic_axes[direction])
        left_d = diffusivity_flat[left_ids]
        right_d = diffusivity_flat[right_ids]
        face_area = np.prod([spacing[i] for i in range(3) if i != direction])
        conductance = _harmonic_mean(left_d, right_d) * face_area / spacing[direction]
        diagonal[left_ids] += conductance
        diagonal[right_ids] += conductance
        row_parts.extend([left_ids, right_ids])
        col_parts.extend([right_ids, left_ids])
        data_parts.extend([-conductance, -conductance])

    # Dirichlet conditions are applied at the external faces of the low/high
    # planes. The half-cell distance produces the correct uniform-medium limit.
    if not morphology.periodic_axes[axis]:
        raise ValueError("the measured axis must use a finite sample with Dirichlet ends")
    # A periodic morphology still represents a finite transport experiment along
    # the measured direction: omit its periodic edge and attach Dirichlet ends.
    ids = np.arange(np.prod(shape), dtype=np.int64).reshape(shape)
    low_slice = [slice(None)] * 3
    high_slice = [slice(None)] * 3
    low_slice[axis] = 0
    high_slice[axis] = shape[axis] - 1
    low_ids = ids[tuple(low_slice)].ravel()
    high_ids = ids[tuple(high_slice)].ravel()
    # The assembly above included the periodic wrap edge. Remove it by rebuilding
    # if the morphology metadata is periodic along the measured axis.
    if morphology.periodic_axes[axis]:
        row_parts = []
        col_parts = []
        data_parts = []
        diagonal[:] = 0.0
        for direction in range(3):
            use_periodic = morphology.periodic_axes[direction] and direction != axis
            left_ids, right_ids = _edge_ids(shape, direction, use_periodic)
            if direction == axis:
                left = [slice(None)] * 3
                right = [slice(None)] * 3
                left[axis] = slice(0, -1)
                right[axis] = slice(1, None)
                left_ids = ids[tuple(left)].ravel()
                right_ids = ids[tuple(right)].ravel()
            left_d = diffusivity_flat[left_ids]
            right_d = diffusivity_flat[right_ids]
            face_area = np.prod([spacing[i] for i in range(3) if i != direction])
            conductance = _harmonic_mean(left_d, right_d) * face_area / spacing[direction]
            diagonal[left_ids] += conductance
            diagonal[right_ids] += conductance
            row_parts.extend([left_ids, right_ids])
            col_parts.extend([right_ids, left_ids])
            data_parts.extend([-conductance, -conductance])

        low_boundary = diffusivity_flat[low_ids] * 2.0 / spacing[axis]
        high_boundary = diffusivity_flat[high_ids] * 2.0 / spacing[axis]
        diagonal[low_ids] += low_boundary
        diagonal[high_ids] += high_boundary
        rhs[low_ids] += low_boundary * low_value
        rhs[high_ids] += high_boundary * high_value
    else:
        # The current implementation requires the measured axis to be finite;
        # this branch supports an explicitly finite morphology for future use.
        raise ValueError("morphology periodicity along measured axis must be True")

    row_parts.append(np.arange(diffusivity.size, dtype=np.int64))
    col_parts.append(np.arange(diffusivity.size, dtype=np.int64))
    data_parts.append(diagonal)
    matrix = coo_matrix(
        (np.concatenate(data_parts), (np.concatenate(row_parts), np.concatenate(col_parts))),
        shape=(diffusivity.size, diffusivity.size),
    ).tocsr()
    return matrix, rhs, diffusivity_flat


def solve_effective_diffusivity(
    morphology: Morphology,
    axis: int,
    d_transport: float = 1.0,
    d_matrix: float = 1.0e-6,
    low_value: float = 1.0,
    high_value: float = 0.0,
    rtol: float = 1.0e-9,
    maxiter: int | None = None,
) -> TransportResult:
    """Solve ``div(D grad c) = 0`` with transverse periodic boundaries.

    The benchmark uses a finite sample in the measured direction with fixed
    concentrations at the two external faces. Interfaces use harmonic face
    conductances, and the matrix coefficient is a documented near-insulating
    regularization rather than an exact zero.
    """

    if axis not in (0, 1, 2):
        raise ValueError("axis must be 0, 1, or 2")
    if d_transport <= 0 or d_matrix <= 0:
        raise ValueError("diffusivities must be positive")
    matrix, rhs, diffusivity = _assemble_system(
        morphology, d_transport, d_matrix, axis, low_value, high_value
    )
    iterations = 0

    def callback(_: np.ndarray) -> None:
        nonlocal iterations
        iterations += 1

    solution, info = cg(matrix, rhs, rtol=rtol, atol=0.0, maxiter=maxiter, callback=callback)
    if info != 0:
        raise RuntimeError(f"conjugate-gradient solver did not converge (info={info})")
    residual = matrix @ solution - rhs
    residual_norm = float(np.linalg.norm(residual) / max(np.linalg.norm(rhs), 1.0))
    concentration = solution.reshape(morphology.shape)
    ids = np.arange(np.prod(morphology.shape), dtype=np.int64).reshape(morphology.shape)
    low_slice = [slice(None)] * 3
    low_slice[axis] = 0
    low_ids = ids[tuple(low_slice)].ravel()
    low_diffusivity = diffusivity[low_ids]
    spacing = morphology.spacing[axis]
    low_fluxes = low_diffusivity * 2.0 / spacing * (low_value - solution[low_ids])
    flux = float(np.mean(low_fluxes))
    length = morphology.physical_size[axis]
    effective = flux * length / (low_value - high_value)
    return TransportResult(axis, float(effective), flux, residual_norm, iterations, concentration)


def analytical_laminate_diffusivity(
    fraction: float, d_transport: float, d_matrix: float, parallel: bool
) -> float:
    """Series/parallel effective diffusivity for a two-phase laminate."""

    if parallel:
        return fraction * d_transport + (1.0 - fraction) * d_matrix
    return 1.0 / (fraction / d_transport + (1.0 - fraction) / d_matrix)
