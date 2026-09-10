"""Benchmark 001: same phase fraction, different morphology and transport."""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from .connectivity import analyze_connectivity, interfacial_area_density, mean_phase_chord_length
from .generators import canonical_morphologies
from .morphology import Morphology
from .transport import FLUX_BALANCE_TOLERANCE, solve_effective_diffusivity
from .validation import run_validation, validation_passes


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _environment() -> dict[str, str]:
    import matplotlib
    import numpy
    import scipy

    return {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": numpy.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
    }


def _source_revision() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL, text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _tracked_source_is_clean() -> bool:
    try:
        status = subprocess.check_output(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            stderr=subprocess.DEVNULL,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return False
    return not status.strip()


def _configuration_digest(configuration: dict[str, Any]) -> str:
    encoded = json.dumps(configuration, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _plot_results(
    morphologies: list[Morphology],
    metric_rows: list[dict[str, Any]],
    transport_rows: list[dict[str, Any]],
    convergence_rows: list[dict[str, Any]],
    contrast_rows: list[dict[str, Any]],
    figure_dir: Path,
) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    by_name = {row["morphology"]: row for row in metric_rows}
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    for column, morphology in enumerate(morphologies[:3]):
        axes[0, column].imshow(morphology.phase[:, :, morphology.shape[2] // 2].T, cmap="viridis")
        axes[0, column].set_title(
            f"{morphology.name}\nf={by_name[morphology.name]['phase_fraction']:.3f}"
        )
        axes[0, column].set_axis_off()
    names = [row["morphology"] for row in metric_rows]
    x = np.arange(len(names))
    width = 0.24
    for offset, axis in enumerate("xyz"):
        values = [
            row["effective_diffusivity"]
            for row in transport_rows
            if row["axis"] == offset
        ]
        axes[1, 0].bar(x + (offset - 1) * width, values, width, label=f"D_eff,{axis}")
    axes[1, 0].set_xticks(x, names, rotation=25, ha="right")
    axes[1, 0].set_ylabel("effective diffusivity")
    axes[1, 0].legend(fontsize=8)
    axes[1, 0].set_title("Directional transport")

    percolation = np.array(
        [[int(row[f"percolates_{axis}"]) for axis in "xyz"] for row in metric_rows]
    )
    axes[1, 1].imshow(percolation, cmap="Greys", vmin=0, vmax=1)
    axes[1, 1].set_xticks(range(3), ["x", "y", "z"])
    axes[1, 1].set_yticks(range(len(names)), names)
    axes[1, 1].set_title("Periodic face-to-face percolation")

    for name in sorted({row["morphology"] for row in convergence_rows}):
        rows = [row for row in convergence_rows if row["morphology"] == name and row["axis"] == 0]
        axes[1, 2].plot(
            [row["resolution"] for row in rows],
            [row["effective_diffusivity"] for row in rows],
            marker="o",
            label=name,
        )
    axes[1, 2].set_xlabel("grid resolution N (N³)")
    axes[1, 2].set_ylabel("D_eff,x")
    axes[1, 2].set_title("Resolution sensitivity")
    axes[1, 2].legend(fontsize=7)
    fig.suptitle("Benchmark 001 — controlled morphology-to-transport comparison")
    fig.savefig(figure_dir / "benchmark_001_overview.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    for name in sorted({row["morphology"] for row in convergence_rows}):
        rows = [row for row in convergence_rows if row["morphology"] == name and row["axis"] == 0]
        ax.plot(
            [row["resolution"] for row in rows],
            [row["effective_diffusivity"] for row in rows],
            marker="o",
            label=name,
        )
    ax.set_xlabel("grid resolution N (N³)")
    ax.set_ylabel("D_eff,x")
    ax.set_title("Benchmark 001 resolution study")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figure_dir / "convergence.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    for name in sorted({row["morphology"] for row in contrast_rows}):
        for direction in range(3):
            rows = [
                row
                for row in contrast_rows
                if row["morphology"] == name and row["axis"] == direction
            ]
            ax.plot(
                [row["d_matrix"] for row in rows],
                [row["effective_diffusivity"] for row in rows],
                marker="o",
                label=f"{name}, {direction}",
            )
    ax.set_xscale("log")
    ax.set_xlabel("matrix diffusivity")
    ax.set_ylabel("effective diffusivity")
    ax.set_title("Transport-contrast sensitivity (N=16)")
    ax.legend(fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(figure_dir / "contrast_sensitivity.png", dpi=180)
    plt.close(fig)


def run_benchmark(
    output_dir: str | Path = "results/benchmark_001",
    resolution: int = 32,
    fraction: float = 0.4,
    d_transport: float = 1.0,
    d_matrix: float = 1.0e-3,
    seed: int = 20260910,
    convergence_resolutions: tuple[int, ...] = (16, 24, 32),
    contrast_values: tuple[float, ...] = (1.0e-2, 1.0e-3, 1.0e-4),
    contrast_resolution: int = 16,
) -> dict[str, Any]:
    """Execute the complete M1 benchmark and write auditable artifacts."""

    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    configuration = {
        "benchmark_id": "benchmark_001",
        "resolution": resolution,
        "fraction": fraction,
        "d_transport": d_transport,
        "d_matrix": d_matrix,
        "seed": seed,
        "convergence_resolutions": list(convergence_resolutions),
        "contrast_values": list(contrast_values),
        "contrast_resolution": contrast_resolution,
        "solver_tolerance": 1.0e-8,
        "flux_balance_tolerance": FLUX_BALANCE_TOLERANCE,
    }
    morphologies = canonical_morphologies((resolution,) * 3, fraction=fraction, seed=seed)
    metric_rows: list[dict[str, Any]] = []
    transport_rows: list[dict[str, Any]] = []
    for morphology in morphologies:
        connectivity = analyze_connectivity(morphology)
        row = {
            "morphology": morphology.name,
            "shape": "x".join(str(value) for value in morphology.shape),
            "target_phase_fraction": fraction,
            "phase_fraction": morphology.volume_fraction,
            "interfacial_area_density": interfacial_area_density(morphology),
            "mean_chord_x": mean_phase_chord_length(morphology, 0),
            "mean_chord_y": mean_phase_chord_length(morphology, 1),
            "mean_chord_z": mean_phase_chord_length(morphology, 2),
            "generator_parameters": json.dumps(morphology.parameters, sort_keys=True),
            **connectivity.as_dict(),
            "content_digest": morphology.content_digest(),
        }
        metric_rows.append(row)
        for axis in range(3):
            result = solve_effective_diffusivity(
                morphology,
                axis=axis,
                d_transport=d_transport,
                d_matrix=d_matrix,
                rtol=1.0e-8,
                maxiter=20000,
            )
            transport_rows.append(
                {
                    "morphology": morphology.name,
                    "axis": axis,
                    "effective_diffusivity": result.effective_diffusivity,
                    "flux": result.flux,
                    "high_flux": result.high_flux,
                    "flux_balance_error": result.flux_balance_error,
                    "flux_balance_tolerance": FLUX_BALANCE_TOLERANCE,
                    "flux_balance_pass": result.flux_balance_error < FLUX_BALANCE_TOLERANCE,
                    "residual_norm": result.residual_norm,
                    "iterations": result.iterations,
                }
            )

    validation = run_validation(shape=(16, 16, 16), d_transport=d_transport, d_matrix=d_matrix)
    convergence_rows: list[dict[str, Any]] = []
    morphology_convergence_rows: list[dict[str, Any]] = []
    for current_resolution in convergence_resolutions:
        current_morphologies = canonical_morphologies(
            (current_resolution,) * 3, fraction=fraction, seed=seed
        )
        for morphology in current_morphologies:
            connectivity = analyze_connectivity(morphology)
            morphology_convergence_rows.append(
                {
                    "morphology": morphology.name,
                    "resolution": current_resolution,
                    "target_phase_fraction": fraction,
                    "phase_fraction": morphology.volume_fraction,
                    "interfacial_area_density": interfacial_area_density(morphology),
                    "mean_chord_x": mean_phase_chord_length(morphology, 0),
                    "mean_chord_y": mean_phase_chord_length(morphology, 1),
                    "mean_chord_z": mean_phase_chord_length(morphology, 2),
                    "generator_parameters": json.dumps(morphology.parameters, sort_keys=True),
                    **connectivity.as_dict(),
                    "content_digest": morphology.content_digest(),
                }
            )
            for axis in range(3):
                result = solve_effective_diffusivity(
                    morphology,
                    axis=axis,
                    d_transport=d_transport,
                    d_matrix=d_matrix,
                    rtol=1.0e-8,
                    maxiter=20000,
                )
                convergence_rows.append(
                    {
                        "morphology": morphology.name,
                        "resolution": current_resolution,
                        "axis": axis,
                        "effective_diffusivity": result.effective_diffusivity,
                        "flux_balance_error": result.flux_balance_error,
                        "residual_norm": result.residual_norm,
                        "iterations": result.iterations,
                    }
                )

    stability_summary: dict[str, Any] = {}
    for morphology_name in {row["morphology"] for row in morphology_convergence_rows}:
        morphology_metrics = [
            row for row in morphology_convergence_rows if row["morphology"] == morphology_name
        ]
        transport_metrics = [
            row for row in convergence_rows if row["morphology"] == morphology_name
        ]
        percolation_states = [
            [row[f"percolates_{axis}"] for axis in "xyz"] for row in morphology_metrics
        ]
        component_counts = [row["components"] for row in morphology_metrics]
        transport_spread = {}
        for axis in range(3):
            values = [
                float(row["effective_diffusivity"])
                for row in transport_metrics
                if row["axis"] == axis
            ]
            transport_spread[f"axis_{axis}"] = (max(values) - min(values)) / max(
                abs(float(np.mean(values))), 1.0e-30
            )
        stability_summary[morphology_name] = {
            "percolation_states_by_resolution": percolation_states,
            "component_counts_by_resolution": component_counts,
            "topology_state_stable": len({json.dumps(state) for state in percolation_states}) == 1
            and len(set(component_counts)) == 1,
            "transport_relative_range_by_axis": transport_spread,
            "phase_fraction_relative_range": (
                max(float(row["phase_fraction"]) for row in morphology_metrics)
                - min(float(row["phase_fraction"]) for row in morphology_metrics)
            )
            / fraction,
            "interfacial_area_density_by_resolution": {
                str(row["resolution"]): row["interfacial_area_density"]
                for row in morphology_metrics
            },
            "mean_chord_lengths_by_resolution": {
                str(row["resolution"]): {
                    axis: row[f"mean_chord_{axis}"] for axis in "xyz"
                }
                for row in morphology_metrics
            },
        }

    contrast_rows: list[dict[str, Any]] = []
    contrast_morphologies = canonical_morphologies(
        (contrast_resolution,) * 3, fraction=fraction, seed=seed
    )
    for current_d_matrix in contrast_values:
        for morphology in contrast_morphologies:
            for axis in range(3):
                result = solve_effective_diffusivity(
                    morphology,
                    axis=axis,
                    d_transport=d_transport,
                    d_matrix=current_d_matrix,
                    rtol=1.0e-8,
                    maxiter=20000,
                )
                contrast_rows.append(
                    {
                        "morphology": morphology.name,
                        "resolution": contrast_resolution,
                        "axis": axis,
                        "d_transport": d_transport,
                        "d_matrix": current_d_matrix,
                        "contrast_ratio": d_transport / current_d_matrix,
                        "effective_diffusivity": result.effective_diffusivity,
                        "flux_balance_error": result.flux_balance_error,
                    }
                )

    _write_csv(root / "morphology_metrics.csv", metric_rows)
    _write_csv(root / "transport_metrics.csv", transport_rows)
    _write_csv(root / "convergence.csv", convergence_rows)
    _write_csv(root / "morphology_convergence.csv", morphology_convergence_rows)
    _write_csv(root / "contrast_sensitivity.csv", contrast_rows)
    (root / "validation.json").write_text(json.dumps(validation, indent=2), encoding="utf-8")
    (root / "stability_summary.json").write_text(
        json.dumps(stability_summary, indent=2), encoding="utf-8"
    )
    manifest = {
        "benchmark_id": "benchmark_001",
        "benchmark_question": (
            "At fixed phase fraction and intrinsic coefficients, "
            "how much can morphology alter transport?"
        ),
        "resolution": resolution,
        "seed": seed,
        "physical_dimensions": [1.0, 1.0, 1.0],
        "voxel_spacing": [1.0 / resolution] * 3,
        "phase_fraction_target": fraction,
        "phase_fraction_actual": {row["morphology"]: row["phase_fraction"] for row in metric_rows},
        "transport_phase": "boolean phase == true",
        "d_transport": d_transport,
        "d_matrix": d_matrix,
        "contrast_sensitivity": {
            "resolution": contrast_resolution,
            "d_matrix_values": list(contrast_values),
            "ratio_values": [d_transport / value for value in contrast_values],
        },
        "boundary_conditions": "Dirichlet c=1/0 along measured axis; periodic transverse axes",
        "solver": "cell-centred finite volume with harmonic face conductance and scipy CG",
        "solver_tolerance": 1.0e-8,
        "flux_balance_tolerance": FLUX_BALANCE_TOLERANCE,
        "validation_reference": "uniform medium and series/parallel laminate limits",
        "runtime_seconds": time.perf_counter() - start,
        "convergence_resolutions": list(convergence_resolutions),
        "configuration": configuration,
        "configuration_digest": _configuration_digest(configuration),
        "environment": _environment(),
        "source_revision": _source_revision(),
        "tracked_source_clean_before_run": _tracked_source_is_clean(),
    }
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (root / "environment.json").write_text(json.dumps(_environment(), indent=2), encoding="utf-8")
    _plot_results(
        morphologies,
        metric_rows,
        transport_rows,
        convergence_rows,
        contrast_rows,
        root / "figures",
    )
    return {
        "manifest": manifest,
        "validation": validation,
        "validation_passed": validation_passes(validation),
        "morphology_metrics": metric_rows,
        "transport_metrics": transport_rows,
        "convergence": convergence_rows,
        "contrast_sensitivity": contrast_rows,
        "stability_summary": stability_summary,
    }
