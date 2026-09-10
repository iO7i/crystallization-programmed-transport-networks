"""Benchmark 002: phenomenological crystalline-like transport perturbations."""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Literal

import matplotlib.pyplot as plt
import numpy as np

from .connectivity import analyze_connectivity, interfacial_area_density
from .crystalline_perturbation import (
    PlacementMode,
    baseline_flux_score,
    build_diffusivity_field,
    select_crystalline_like,
)
from .generators import canonical_morphologies
from .morphology import Morphology
from .transport import (
    FLUX_BALANCE_TOLERANCE,
    solve_effective_diffusivity,
    solve_effective_diffusivity_field,
)

B002_SOLVER_RTOL = 1.0e-10
B002_RELATIVE_TOLERANCE = 1.0e-10
B002_ABSOLUTE_TOLERANCE = 1.0e-12
STABLE_MORPHOLOGY_NAMES = ("lamellae", "cylinders", "gyroid_level_set")
DEFAULT_FRACTIONS = (0.0, 0.1, 0.2, 0.4, 0.6)
DEFAULT_MOBILITY_RATIOS = (1.0, 0.1, 0.01, 0.001)
DEFAULT_SWEEP_SEEDS = (20260910, 20260911, 20260912)
DEFAULT_PLACEMENT_SEEDS = (20260910, 20260911, 20260912, 20260913, 20260914)
PLACEMENT_MODES: tuple[PlacementMode, ...] = (
    "random",
    "interface",
    "baseline-flux-ranked",
    "low_criticality",
)
ValidationProfile = Literal["canonical", "smoke"]


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


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


def _stable_morphologies(resolution: int, fraction: float, seed: int) -> list[Morphology]:
    return canonical_morphologies(
        (resolution, resolution, resolution), fraction=fraction, seed=seed
    )[:3]


def _baseline(
    morphology: Morphology,
    d_mobile: float,
    d_matrix: float,
) -> tuple[Any, Any, Any]:
    return tuple(
        solve_effective_diffusivity(
            morphology,
            axis=axis,
            d_transport=d_mobile,
            d_matrix=d_matrix,
            rtol=B002_SOLVER_RTOL,
            maxiter=20000,
        )
        for axis in range(3)
    )  # type: ignore[return-value]


def _cache_entry(
    morphology: Morphology,
    d_mobile: float,
    d_matrix: float,
) -> tuple[tuple[Any, Any, Any], np.ndarray]:
    baseline = _baseline(morphology, d_mobile, d_matrix)
    score = baseline_flux_score(morphology, baseline, d_mobile, d_matrix)
    return baseline, score


def _condition_row(
    morphology: Morphology,
    baseline: tuple[Any, Any, Any],
    baseline_flux: np.ndarray,
    fraction_of_transport_phase: float,
    placement_mode: PlacementMode,
    seed: int,
    d_mobile: float,
    d_crystal: float,
    d_matrix: float,
    experiment: str,
) -> dict[str, Any]:
    placement = select_crystalline_like(
        morphology,
        fraction_of_transport_phase,
        placement_mode,
        seed,
        baseline_flux=baseline_flux,
    )
    field = build_diffusivity_field(
        morphology,
        placement.crystalline_mask,
        d_mobile=d_mobile,
        d_crystal=d_crystal,
        d_matrix=d_matrix,
    )
    results = [
        solve_effective_diffusivity_field(
            morphology,
            field,
            axis=axis,
            rtol=B002_SOLVER_RTOL,
            maxiter=20000,
        )
        for axis in range(3)
    ]
    baseline_values = [float(result.effective_diffusivity) for result in baseline]
    effective_values = [float(result.effective_diffusivity) for result in results]
    normalized_values = [
        value / baseline_value if baseline_value else float("nan")
        for value, baseline_value in zip(effective_values, baseline_values, strict=True)
    ]
    connectivity = analyze_connectivity(morphology).as_dict()
    row: dict[str, Any] = {
        "experiment": experiment,
        "morphology": morphology.name,
        "resolution": morphology.shape[0],
        "transport_phase_fraction": morphology.volume_fraction,
        "crystalline_like_fraction_definition": "fraction of transport phase voxels",
        "crystalline_like_fraction_target": fraction_of_transport_phase,
        "crystalline_like_fraction_actual_transport_phase": (
            placement.actual_fraction_of_transport_phase
        ),
        "crystalline_like_fraction_actual_total_material": (
            placement.actual_total_material_fraction
        ),
        "placement_mode": placement.mode,
        "seed": seed,
        "eligible_voxel_count": placement.eligible_voxel_count,
        "selected_voxel_count": placement.selected_voxel_count,
        "placement_algorithm": placement.algorithm,
        "selection_digest": placement.selection_digest,
        "d_mobile": d_mobile,
        "d_crystal": d_crystal,
        "d_matrix": d_matrix,
        "contrast_ratio_mobile_to_matrix": d_mobile / d_matrix,
        "geometry_perturbation_mode": "none",
        **connectivity,
    }
    for axis, label in enumerate("xyz"):
        row[f"d_eff_{label}"] = effective_values[axis]
        row[f"d_eff_normalized_{label}"] = normalized_values[axis]
        row[f"d_eff_baseline_{label}"] = baseline_values[axis]
        row[f"flux_error_{label}"] = results[axis].flux_balance_error
        row[f"residual_norm_{label}"] = results[axis].residual_norm
        row[f"iterations_{label}"] = results[axis].iterations
    return row


def _aggregate_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not rows:
        return []
    keys = (
        "experiment",
        "morphology",
        "resolution",
        "crystalline_like_fraction_target",
        "placement_mode",
        "d_mobile",
        "d_crystal",
        "d_matrix",
    )
    groups: dict[tuple[object, ...], list[dict[str, Any]]] = {}
    for row in rows:
        key = tuple(row[key_name] for key_name in keys)
        groups.setdefault(key, []).append(row)
    output: list[dict[str, Any]] = []
    for key, group in sorted(groups.items(), key=lambda item: tuple(str(v) for v in item[0])):
        row = {key_name: value for key_name, value in zip(keys, key, strict=True)}
        row["seed_count"] = len(group)
        row["crystalline_like_fraction_actual_transport_phase_mean"] = float(
            np.mean(
                [float(item["crystalline_like_fraction_actual_transport_phase"]) for item in group]
            )
        )
        row["crystalline_like_fraction_actual_transport_phase_std"] = float(
            np.std(
                [float(item["crystalline_like_fraction_actual_transport_phase"]) for item in group]
            )
        )
        for label in "xyz":
            values = np.asarray([float(item[f"d_eff_{label}"]) for item in group])
            normalized = np.asarray([float(item[f"d_eff_normalized_{label}"]) for item in group])
            row[f"d_eff_mean_{label}"] = float(np.mean(values))
            row[f"d_eff_std_{label}"] = float(np.std(values))
            row[f"d_eff_normalized_mean_{label}"] = float(np.mean(normalized))
            row[f"d_eff_normalized_std_{label}"] = float(np.std(normalized))
            row[f"max_flux_error_{label}"] = max(
                float(item[f"flux_error_{label}"]) for item in group
            )
        row["percolates_x"] = group[0]["percolates_x"]
        row["percolates_y"] = group[0]["percolates_y"]
        row["percolates_z"] = group[0]["percolates_z"]
        row["components"] = group[0]["components"]
        row["largest_component_fraction"] = group[0]["largest_component_fraction"]
        output.append(row)
    return output


def _max_relative_error(reference: float, value: float) -> float:
    return abs(value - reference) / max(abs(reference), 1.0e-30)


def _load_benchmark001_reference(reference_dir: Path) -> tuple[str, list[dict[str, str]]]:
    manifest = json.loads((reference_dir / "manifest.json").read_text(encoding="utf-8"))
    rows = _read_csv(reference_dir / "transport_metrics.csv")
    return str(manifest["source_revision"]), rows


def _control_checks(
    root: Path,
    cache: dict[tuple[int, str, float], tuple[tuple[Any, Any, Any], np.ndarray, Morphology]],
    sweep_rows: list[dict[str, Any]],
    reference_dir: Path,
    d_mobile: float,
    d_matrix: float,
    fraction: float,
    seed: int,
) -> dict[str, Any]:
    reference_source_revision, reference_rows = _load_benchmark001_reference(reference_dir)
    canonical_errors: list[float] = []
    for key, (baseline, _, _morphology) in cache.items():
        resolution, morphology_name, current_d_matrix = key
        if resolution != 32 or current_d_matrix != d_matrix:
            continue
        reference_for_morphology = [
            row for row in reference_rows if row["morphology"] == morphology_name
        ]
        for axis, result in enumerate(baseline):
            reference = next(row for row in reference_for_morphology if int(row["axis"]) == axis)
            canonical_errors.append(
                _max_relative_error(
                    float(reference["effective_diffusivity"]),
                    float(result.effective_diffusivity),
                )
            )
    no_op_rows = [
        row
        for row in sweep_rows
        if float(row["d_crystal"]) == d_mobile
    ]
    no_op_errors = [
        abs(float(row[f"d_eff_normalized_{axis}"]) - 1.0)
        for row in no_op_rows
        for axis in "xyz"
    ]
    sample_key = next(
        key
        for key in sorted(cache)
        if key[1] == "gyroid_level_set" and key[2] == d_matrix
    )
    sample_baseline, sample_flux, sample_morphology = cache[sample_key]
    first = select_crystalline_like(sample_morphology, 0.2, "random", seed, sample_flux)
    second = select_crystalline_like(sample_morphology, 0.2, "random", seed, sample_flux)
    label_field = build_diffusivity_field(
        sample_morphology,
        first.crystalline_mask,
        d_mobile,
        d_mobile,
        d_matrix,
    )
    baseline_field = np.where(sample_morphology.phase, d_mobile, d_matrix).astype(float)
    field_solver_errors = []
    for axis, scalar_result in enumerate(sample_baseline):
        field_result = solve_effective_diffusivity_field(
            sample_morphology,
            baseline_field,
            axis,
            rtol=B002_SOLVER_RTOL,
            maxiter=20000,
        )
        field_solver_errors.append(
            _max_relative_error(
                scalar_result.effective_diffusivity,
                field_result.effective_diffusivity,
            )
        )
    return {
        "benchmark001_reference_source_revision": reference_source_revision,
        "benchmark001_baseline_max_relative_error": max(canonical_errors, default=float("inf")),
        "benchmark001_baseline_tolerance": B002_RELATIVE_TOLERANCE,
        "benchmark001_baseline_pass": max(canonical_errors, default=float("inf"))
        <= B002_RELATIVE_TOLERANCE,
        "d_crystal_equals_d_mobile_noop_max_absolute_normalized_error": max(
            no_op_errors, default=float("inf")
        ),
        "d_crystal_equals_d_mobile_noop_tolerance": B002_RELATIVE_TOLERANCE,
        "d_crystal_equals_d_mobile_noop_pass": max(no_op_errors, default=float("inf"))
        <= B002_RELATIVE_TOLERANCE,
        "same_seed_selection_digest_equal": first.selection_digest == second.selection_digest,
        "label_only_noop_field_exact": bool(np.array_equal(label_field, baseline_field)),
        "field_solver_baseline_max_relative_error": max(field_solver_errors, default=float("inf")),
        "field_solver_baseline_pass": max(field_solver_errors, default=float("inf"))
        <= B002_RELATIVE_TOLERANCE,
        "reference_result_directory": str(root / ".." / reference_dir),
        "zero_fraction_definition": "zero selected transport-phase voxels",
        "matrix_diffusivity_explicit": d_matrix,
        "fraction_parameter_used_for_sample": fraction,
    }


def _smoke_checks(
    root: Path,
    sample_morphology: Morphology,
    masks: dict[str, np.ndarray],
    sweep_rows: list[dict[str, Any]],
    placement_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    """Run fast, exact/core checks for an explicitly labelled CI smoke profile."""

    required_outputs = (
        "runs.csv",
        "summary.csv",
        "placement_summary.csv",
        "grid_summary.csv",
        "contrast_sensitivity.csv",
        "morphology_metrics.csv",
        "placement_masks.npz",
        "figures/mobility_sweep.png",
        "figures/placement_comparison.png",
        "figures/placement_slices.png",
    )
    artifact_generation_pass = all(
        (root / relative_path).is_file() for relative_path in required_outputs
    )

    rows = sweep_rows + placement_rows
    denominator_checks = []
    for row in rows:
        eligible = int(row["eligible_voxel_count"])
        selected = int(row["selected_voxel_count"])
        target = float(row["crystalline_like_fraction_target"])
        actual = float(row["crystalline_like_fraction_actual_transport_phase"])
        denominator_checks.append(
            eligible > 0
            and 0 <= selected <= eligible
            and selected == round(target * eligible)
            and np.isclose(actual, selected / eligible, rtol=1.0e-15, atol=1.0e-15)
        )
    phase_denominator_pass = bool(denominator_checks) and all(denominator_checks)

    expected_modes = set(PLACEMENT_MODES)
    valid_masks_pass = (
        set(masks) == {"phase", *expected_modes}
        and masks["phase"].shape == sample_morphology.shape
        and masks["phase"].dtype == bool
        and all(
            mask.shape == sample_morphology.shape
            and mask.dtype == bool
            and np.all(~mask | sample_morphology.phase)
            for name, mask in masks.items()
            if name != "phase"
        )
    )

    invalid_fraction_rejected = False
    try:
        select_crystalline_like(sample_morphology, 1.1, "random", 0)
    except ValueError:
        invalid_fraction_rejected = True

    invalid_mask_rejected = False
    matrix_indices = np.flatnonzero(~sample_morphology.phase.ravel())
    if matrix_indices.size:
        invalid_mask = np.zeros(sample_morphology.shape, dtype=bool)
        invalid_mask.ravel()[matrix_indices[0]] = True
        try:
            build_diffusivity_field(sample_morphology, invalid_mask, 1.0, 0.1, 1.0e-3)
        except ValueError:
            invalid_mask_rejected = True

    invalid_field_rejected = False
    try:
        solve_effective_diffusivity_field(
            sample_morphology,
            np.zeros(sample_morphology.shape, dtype=float),
            axis=0,
        )
    except ValueError:
        invalid_field_rejected = True

    fail_closed_inputs_pass = (
        invalid_fraction_rejected and invalid_mask_rejected and invalid_field_rejected
    )
    return {
        "artifact_generation_pass": artifact_generation_pass,
        "phase_denominator_pass": phase_denominator_pass,
        "valid_masks_pass": valid_masks_pass,
        "fail_closed_inputs_pass": fail_closed_inputs_pass,
    }


def _grid_stability(grid_rows: list[dict[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {
        "resolutions": sorted({int(row["resolution"]) for row in grid_rows}),
        "cases": {},
        "structural_connectivity_transition": {
            "tested": False,
            "reason": "002C structural perturbations are deferred until 002A and 002B are stable",
        },
    }
    for case in sorted({str(row["experiment"]) for row in grid_rows}):
        case_rows = [row for row in grid_rows if row["experiment"] == case]
        case_result: dict[str, Any] = {}
        for morphology in sorted({str(row["morphology"]) for row in case_rows}):
            morphology_rows = [row for row in case_rows if row["morphology"] == morphology]
            mode_groups = (
                sorted({str(row["placement_mode"]) for row in morphology_rows})
                if case == "grid_placement_comparison"
                else ["all"]
            )
            morphology_result: dict[str, Any] = {}
            for mode in mode_groups:
                selected_rows = (
                    morphology_rows
                    if mode == "all"
                    else [row for row in morphology_rows if row["placement_mode"] == mode]
                )
                topology = [
                    [row[f"percolates_{axis}"] for axis in "xyz"] + [row["components"]]
                    for row in selected_rows
                ]
                spreads: dict[str, float] = {}
                for axis in "xyz":
                    values = [float(row[f"d_eff_{axis}"]) for row in selected_rows]
                    spreads[axis] = (max(values) - min(values)) / max(
                        abs(float(np.mean(values))), 1.0e-30
                    )
                morphology_result[mode] = {
                    "topology_states_by_resolution": topology,
                    "topology_state_stable": len({json.dumps(item) for item in topology}) == 1,
                    "effective_transport_relative_range": spreads,
                }
            case_result[morphology] = morphology_result
        result["cases"][case] = case_result
    return result


def _plot_figures(root: Path) -> None:
    summary_rows = _read_csv(root / "summary.csv")
    placement_rows = _read_csv(root / "placement_summary.csv")
    fig_dir = root / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)
    selected_ratio = 0.1
    for axis, label in enumerate("xyz"):
        for morphology in STABLE_MORPHOLOGY_NAMES:
            rows = [
                row
                for row in summary_rows
                if row["morphology"] == morphology
                and abs(float(row["d_crystal"]) - selected_ratio) < 1e-15
            ]
            rows.sort(key=lambda row: float(row["crystalline_like_fraction_target"]))
            axes[axis].errorbar(
                [float(row["crystalline_like_fraction_target"]) for row in rows],
                [float(row[f"d_eff_normalized_mean_{label}"]) for row in rows],
                yerr=[float(row[f"d_eff_normalized_std_{label}"]) for row in rows],
                marker="o",
                capsize=2,
                label=morphology,
            )
        axes[axis].set_xlabel("crystalline-like fraction of transport phase")
        axes[axis].set_ylabel(f"normalized D_eff,{label}")
        axes[axis].set_ylim(bottom=0)
        axes[axis].set_title(f"D_crystal / D_mobile = {selected_ratio:g}")
    axes[0].legend(fontsize=7)
    fig.suptitle("Benchmark 002A — mobility-only perturbation")
    fig.savefig(fig_dir / "mobility_sweep.png", dpi=180)
    plt.close(fig)

    modes = [mode for mode in PLACEMENT_MODES]
    x = np.arange(len(modes))
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    for axis, label in enumerate("xyz"):
        rows = [row for row in placement_rows if row["placement_mode"] in modes]
        values = {
            row["placement_mode"]: float(row[f"d_eff_normalized_mean_{label}"])
            for row in rows
        }
        errors = {
            row["placement_mode"]: float(row[f"d_eff_normalized_std_{label}"])
            for row in rows
        }
        axes[axis].bar(
            x,
            [values[mode] for mode in modes],
            yerr=[errors[mode] for mode in modes],
            capsize=3,
        )
        axes[axis].set_xticks(x, modes, rotation=25, ha="right")
        axes[axis].set_ylabel(f"normalized D_eff,{label}")
        axes[axis].set_ylim(bottom=0)
    fig.suptitle("Benchmark 002B — equal fraction, different placement (gyroid)")
    fig.savefig(fig_dir / "placement_comparison.png", dpi=180)
    plt.close(fig)

    with np.load(root / "placement_masks.npz", allow_pickle=False) as masks:
        phase = masks["phase"]
        names = [
            ("random", "random"),
            ("interface", "interface"),
            ("baseline_flux_ranked", "baseline-flux-ranked"),
            ("low_criticality", "low_criticality"),
        ]
        fig, axes = plt.subplots(1, 4, figsize=(12, 3), constrained_layout=True)
        for axis, (archive_name, display_name) in enumerate(names):
            mask = masks[archive_name].astype(bool)
            rgb = np.zeros((*phase.shape[:2], 3), dtype=float)
            slice_phase = phase[:, :, phase.shape[2] // 2]
            slice_mask = mask[:, :, mask.shape[2] // 2]
            rgb[slice_phase] = (0.18, 0.65, 0.25)
            rgb[~slice_phase] = (0.05, 0.05, 0.05)
            rgb[slice_mask] = (0.95, 0.75, 0.10)
            axes[axis].imshow(rgb.transpose(1, 0, 2), origin="lower")
            axes[axis].set_title(display_name)
            axes[axis].set_axis_off()
        fig.suptitle("Benchmark 002B — gyroid placement sanity slices")
        fig.savefig(fig_dir / "placement_slices.png", dpi=180)
        plt.close(fig)


def run_benchmark(
    output_dir: str | Path = "results/benchmark_002",
    resolution: int = 16,
    fraction: float = 0.4,
    d_mobile: float = 1.0,
    d_matrix: float = 1.0e-3,
    crystalline_fractions: tuple[float, ...] = DEFAULT_FRACTIONS,
    mobility_ratios: tuple[float, ...] = DEFAULT_MOBILITY_RATIOS,
    sweep_seeds: tuple[int, ...] = DEFAULT_SWEEP_SEEDS,
    placement_seeds: tuple[int, ...] = DEFAULT_PLACEMENT_SEEDS,
    grid_resolutions: tuple[int, ...] = (16, 24, 32),
    contrast_matrix_values: tuple[float, ...] = (1.0e-2, 1.0e-3, 1.0e-4),
    benchmark001_reference_dir: str | Path = "results/benchmark_001",
    validation_profile: ValidationProfile = "canonical",
) -> dict[str, Any]:
    """Run bounded 002A/002B experiments; structural 002C is intentionally deferred."""

    if validation_profile not in {"canonical", "smoke"}:
        raise ValueError("validation_profile must be 'canonical' or 'smoke'")

    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    reference_dir = Path(benchmark001_reference_dir)
    start = time.perf_counter()
    configuration = {
        "benchmark_id": "benchmark_002",
        "source_benchmark": "benchmark_001",
        "resolution": resolution,
        "transport_phase_fraction": fraction,
        "d_mobile": d_mobile,
        "d_crystal_ratios": list(mobility_ratios),
        "d_matrix": d_matrix,
        "crystalline_fractions_of_transport_phase": list(crystalline_fractions),
        "sweep_seeds": list(sweep_seeds),
        "placement_seeds": list(placement_seeds),
        "placement_modes": list(PLACEMENT_MODES),
        "grid_resolutions": list(grid_resolutions),
        "contrast_matrix_values": list(contrast_matrix_values),
        "solver_tolerance": B002_SOLVER_RTOL,
        "flux_balance_tolerance": FLUX_BALANCE_TOLERANCE,
        "geometry_perturbation_mode": "none",
        "validation_profile": validation_profile,
    }
    cache: dict[tuple[int, str, float], tuple[tuple[Any, Any, Any], np.ndarray, Morphology]] = {}

    def get_cache(
        current_resolution: int, name: str, current_d_matrix: float
    ) -> tuple[Any, Any, Morphology]:
        key = (current_resolution, name, current_d_matrix)
        if key not in cache:
            morphologies = _stable_morphologies(current_resolution, fraction, sweep_seeds[0])
            morphology = next(item for item in morphologies if item.name == name)
            baseline, flux = _cache_entry(morphology, d_mobile, current_d_matrix)
            cache[key] = (baseline, flux, morphology)
        return cache[key]

    sweep_rows: list[dict[str, Any]] = []
    for morphology in _stable_morphologies(resolution, fraction, sweep_seeds[0]):
        baseline, baseline_flux, _ = get_cache(resolution, morphology.name, d_matrix)
        for current_fraction in crystalline_fractions:
            for ratio in mobility_ratios:
                for seed in sweep_seeds:
                    sweep_rows.append(
                        _condition_row(
                            morphology,
                            baseline,
                            baseline_flux,
                            current_fraction,
                            "random",
                            seed,
                            d_mobile,
                            d_mobile * ratio,
                            d_matrix,
                            "002A_mobility_sweep",
                        )
                    )

    placement_rows: list[dict[str, Any]] = []
    gyroid = next(
        morphology
        for morphology in _stable_morphologies(resolution, fraction, sweep_seeds[0])
        if morphology.name == "gyroid_level_set"
    )
    gyroid_baseline, gyroid_flux, _ = get_cache(resolution, gyroid.name, d_matrix)
    for mode in PLACEMENT_MODES:
        for seed in placement_seeds:
            placement_rows.append(
                _condition_row(
                    gyroid,
                    gyroid_baseline,
                    gyroid_flux,
                    0.2,
                    mode,
                    seed,
                    d_mobile,
                    d_mobile * 0.1,
                    d_matrix,
                    "002B_placement",
                )
            )

    contrast_rows: list[dict[str, Any]] = []
    for current_d_matrix in contrast_matrix_values:
        baseline, baseline_flux, morphology = get_cache(
            resolution, "gyroid_level_set", current_d_matrix
        )
        for mode in PLACEMENT_MODES:
            contrast_rows.append(
                _condition_row(
                    morphology,
                    baseline,
                    baseline_flux,
                    0.2,
                    mode,
                    placement_seeds[0],
                    d_mobile,
                    d_mobile * 0.1,
                    current_d_matrix,
                    "002B_contrast",
                )
            )

    grid_rows: list[dict[str, Any]] = []
    for current_resolution in grid_resolutions:
        for morphology in _stable_morphologies(current_resolution, fraction, sweep_seeds[0]):
            baseline, baseline_flux, _ = get_cache(current_resolution, morphology.name, d_matrix)
            grid_rows.append(
                _condition_row(
                    morphology,
                    baseline,
                    baseline_flux,
                    0.0,
                    "random",
                    sweep_seeds[0],
                    d_mobile,
                    d_mobile * 0.1,
                    d_matrix,
                    "grid_baseline",
                )
            )
            grid_rows.append(
                _condition_row(
                    morphology,
                    baseline,
                    baseline_flux,
                    0.2,
                    "random",
                    sweep_seeds[0],
                    d_mobile,
                    d_mobile * 0.1,
                    d_matrix,
                    "grid_moderate_random",
                )
            )
        gyroid_grid = next(
            item
            for item in _stable_morphologies(current_resolution, fraction, sweep_seeds[0])
            if item.name == "gyroid_level_set"
        )
        baseline, baseline_flux, _ = get_cache(current_resolution, gyroid_grid.name, d_matrix)
        for mode in PLACEMENT_MODES:
            grid_rows.append(
                _condition_row(
                    gyroid_grid,
                    baseline,
                    baseline_flux,
                    0.2,
                    mode,
                    sweep_seeds[0],
                    d_mobile,
                    d_mobile * 0.1,
                    d_matrix,
                    "grid_placement_comparison",
                )
            )

    morphology_rows: list[dict[str, Any]] = []
    for morphology in _stable_morphologies(resolution, fraction, sweep_seeds[0]):
        morphology_rows.append(
            {
                "morphology": morphology.name,
                "resolution": resolution,
                "transport_phase_fraction": morphology.volume_fraction,
                "interfacial_area_density": interfacial_area_density(morphology),
                "generator_parameters": json.dumps(morphology.parameters, sort_keys=True),
                **analyze_connectivity(morphology).as_dict(),
            }
        )

    all_rows = sweep_rows + placement_rows + contrast_rows + grid_rows
    summary_rows = _aggregate_rows(sweep_rows)
    placement_summary_rows = _aggregate_rows(placement_rows)
    grid_summary_rows = _aggregate_rows(grid_rows)
    _write_csv(root / "runs.csv", all_rows)
    _write_csv(root / "summary.csv", summary_rows)
    _write_csv(root / "placement_summary.csv", placement_summary_rows)
    _write_csv(root / "grid_summary.csv", grid_summary_rows)
    _write_csv(root / "contrast_sensitivity.csv", contrast_rows)
    _write_csv(root / "morphology_metrics.csv", morphology_rows)

    mask_source = next(
        item
        for item in _stable_morphologies(resolution, fraction, sweep_seeds[0])
        if item.name == "gyroid_level_set"
    )
    mask_baseline, mask_flux, _ = get_cache(resolution, mask_source.name, d_matrix)
    masks: dict[str, np.ndarray] = {"phase": mask_source.phase}
    for mode in PLACEMENT_MODES:
        masks[mode] = select_crystalline_like(
            mask_source,
            0.2,
            mode,
            placement_seeds[0],
            baseline_flux=mask_flux,
        ).crystalline_mask
    np.savez_compressed(
        root / "placement_masks.npz",
        phase=masks["phase"],
        random=masks["random"],
        interface=masks["interface"],
        baseline_flux_ranked=masks["baseline-flux-ranked"],
        low_criticality=masks["low_criticality"],
    )

    for morphology_name in STABLE_MORPHOLOGY_NAMES:
        get_cache(32, morphology_name, d_matrix)

    controls = _control_checks(
        root,
        cache,
        sweep_rows,
        reference_dir,
        d_mobile,
        d_matrix,
        fraction,
        sweep_seeds[0],
    )
    _plot_figures(root)
    smoke_controls = _smoke_checks(root, mask_source, masks, sweep_rows, placement_rows)
    smoke_controls.update(
        {
            "deterministic_placement_pass": bool(
                controls["same_seed_selection_digest_equal"]
            ),
            "no_op_equivalence_pass": bool(
                controls["d_crystal_equals_d_mobile_noop_pass"]
                and controls["label_only_noop_field_exact"]
            ),
        }
    )
    max_flux_error = max(
        float(row[f"flux_error_{axis}"])
        for row in all_rows
        for axis in "xyz"
    )
    canonical_numerical_acceptance = bool(
        all(
            bool(value)
            for key, value in controls.items()
            if key.endswith("_pass") or key.endswith("_equal") or key.endswith("_exact")
        )
        and max_flux_error < FLUX_BALANCE_TOLERANCE
    )
    smoke_core_acceptance = bool(all(bool(value) for value in smoke_controls.values()))
    validation = {
        "validation_profile": validation_profile,
        "validation_passed": (
            canonical_numerical_acceptance and smoke_core_acceptance
            if validation_profile == "canonical"
            else smoke_core_acceptance
        ),
        "flux_balance_tolerance": FLUX_BALANCE_TOLERANCE,
        "maximum_flux_balance_error": max_flux_error,
        "flux_balance_acceptance_applied": validation_profile == "canonical",
        "canonical_numerical_acceptance": canonical_numerical_acceptance,
        "smoke_core_acceptance": smoke_core_acceptance,
        "controls": controls,
        "smoke_controls": smoke_controls,
        "stochastic_replicates": {
            "mobility_sweep_seed_count": len(sweep_seeds),
            "placement_comparison_seed_count": len(placement_seeds),
            "placement_seed_policy": "fixed explicit seeds recorded per run",
        },
    }
    (root / "validation.json").write_text(json.dumps(validation, indent=2), encoding="utf-8")
    stability = {
        "grid_stability": _grid_stability(grid_rows),
        "placement_summary_file": "placement_summary.csv",
        "contrast_summary_file": "contrast_sensitivity.csv",
        "random_field_excluded": True,
        "structural_perturbation_status": "deferred",
    }
    (root / "stability_summary.json").write_text(json.dumps(stability, indent=2), encoding="utf-8")
    b001_manifest = json.loads((reference_dir / "manifest.json").read_text(encoding="utf-8"))
    b001_audit_path = reference_dir / "reproduction_audit.json"
    b001_audit = (
        json.loads(b001_audit_path.read_text(encoding="utf-8")) if b001_audit_path.exists() else {}
    )
    manifest = {
        "benchmark_id": "benchmark_002",
        "title": "When Crystallization Helps or Hurts Transport",
        "scientific_boundary": "phenomenological crystallinity-perturbation benchmark",
        "question": (
            "How do controlled low-mobility labels and their spatial placement affect "
            "transport in already-defined synthetic networks?"
        ),
        "source_benchmark": "benchmark_001",
        "starting_morphologies": list(STABLE_MORPHOLOGY_NAMES),
        "excluded_morphologies": ["gaussian_random_field"],
        "transport_phase_fraction": fraction,
        "crystalline_fraction_definition": (
            "fraction of existing transport-phase voxels labelled crystalline-like"
        ),
        "d_mobile": d_mobile,
        "d_crystal_ratios": list(mobility_ratios),
        "d_matrix": d_matrix,
        "d_matrix_contrast_values": list(contrast_matrix_values),
        "placement_modes": list(PLACEMENT_MODES),
        "geometry_perturbation_mode": "none; 002C deferred",
        "resolution": resolution,
        "grid_resolutions": list(grid_resolutions),
        "solver": "Benchmark 001 cell-centred finite volume with voxelwise D(x) field",
        "solver_tolerance": B002_SOLVER_RTOL,
        "flux_balance_tolerance": FLUX_BALANCE_TOLERANCE,
        "boundary_conditions": "Dirichlet c=1/0 along measured axis; periodic transverse axes",
        "runtime_seconds": time.perf_counter() - start,
        "configuration": configuration,
        "configuration_digest": _configuration_digest(configuration),
        "environment": _environment(),
        "source_revision": _source_revision(),
        "tracked_source_clean_before_run": _tracked_source_is_clean(),
        "benchmark001_reference_source_revision": b001_manifest["source_revision"],
        "benchmark001_reference_result_revision": b001_audit.get(
            "reference_result_commit", "unknown"
        ),
        "benchmark001_reference_manifest": str(reference_dir / "manifest.json"),
        "headline_interpretation": (
            "Fixed-geometry passive mobility perturbation control; 002C structural "
            "help/hurt remains open and deferred"
        ),
        "validation_profile": validation_profile,
    }
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (root / "environment.json").write_text(json.dumps(_environment(), indent=2), encoding="utf-8")
    return {
        "validation_passed": validation["validation_passed"],
        "validation": validation,
        "validation_profile": validation_profile,
        "manifest": manifest,
        "summary": summary_rows,
        "placement_summary": placement_summary_rows,
        "stability": stability,
    }
