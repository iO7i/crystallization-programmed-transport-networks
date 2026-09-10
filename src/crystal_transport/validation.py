"""Analytical validation experiments for the diffusion solver."""

from __future__ import annotations

from typing import Any, cast

from .generators import lamellae
from .morphology import Morphology
from .transport import (
    FLUX_BALANCE_TOLERANCE,
    analytical_laminate_diffusivity,
    solve_effective_diffusivity,
)


def _uniform(shape: tuple[int, int, int], value: bool = True) -> Morphology:
    import numpy as np

    return Morphology(np.full(shape, value, dtype=bool), name="uniform")


def run_validation(
    shape: tuple[int, int, int] = (16, 16, 16),
    fraction: float = 0.4,
    d_transport: float = 1.0,
    d_matrix: float = 1.0e-3,
) -> dict[str, object]:
    """Run uniform and laminate-limit checks without hard-coding solver output."""

    uniform = _uniform(shape)
    uniform_results: list[dict[str, Any]] = []
    for axis in range(3):
        result = solve_effective_diffusivity(
            uniform, axis=axis, d_transport=d_transport, d_matrix=d_transport
        )
        uniform_results.append(
            {
                "axis": axis,
                "numerical": result.effective_diffusivity,
                "expected": d_transport,
                "flux_balance_error": result.flux_balance_error,
            }
        )

    # Use an exactly representable 50/50 laminate for the analytical limit;
    # Benchmark 001 separately reports the phase-fraction discretization it uses.
    laminate = lamellae(shape, fraction=0.5, axis=0, period_voxels=2)
    transverse = solve_effective_diffusivity(
        laminate, axis=1, d_transport=d_transport, d_matrix=d_matrix
    )
    normal = solve_effective_diffusivity(
        laminate, axis=0, d_transport=d_transport, d_matrix=d_matrix
    )
    expected_parallel = analytical_laminate_diffusivity(
        laminate.volume_fraction, d_transport, d_matrix, parallel=True
    )
    expected_series = analytical_laminate_diffusivity(
        laminate.volume_fraction, d_transport, d_matrix, parallel=False
    )
    validation_result: dict[str, Any] = {
        "uniform": uniform_results,
        "laminate": [
            {
                "axis": 0,
                "numerical": normal.effective_diffusivity,
                "expected": expected_series,
                "relative_error": abs(normal.effective_diffusivity - expected_series)
                / expected_series,
                "flux_balance_error": normal.flux_balance_error,
            },
            {
                "axis": 1,
                "numerical": transverse.effective_diffusivity,
                "expected": expected_parallel,
                "relative_error": abs(transverse.effective_diffusivity - expected_parallel)
                / expected_parallel,
                "flux_balance_error": transverse.flux_balance_error,
            },
        ],
        "parameters": {
            "shape": list(shape),
            "fraction_requested": fraction,
            "validation_laminate_fraction": laminate.volume_fraction,
            "fraction_actual": laminate.volume_fraction,
            "d_transport": d_transport,
            "d_matrix": d_matrix,
        },
    }
    return validation_result


def validation_passes(result: dict[str, object], tolerance: float = 2.0e-2) -> bool:
    uniform = cast(list[dict[str, Any]], result["uniform"])
    laminate = cast(list[dict[str, Any]], result["laminate"])
    uniform_ok = all(
        abs(float(item["numerical"]) - float(item["expected"])) / float(item["expected"])
        < tolerance
        for item in uniform
    )
    uniform_flux_ok = all(
        float(item["flux_balance_error"]) < FLUX_BALANCE_TOLERANCE
        for item in uniform
    )
    laminate_ok = all(
        float(item["relative_error"]) < tolerance
        and float(item["flux_balance_error"]) < FLUX_BALANCE_TOLERANCE
        for item in laminate
    )
    return uniform_ok and uniform_flux_ok and laminate_ok
