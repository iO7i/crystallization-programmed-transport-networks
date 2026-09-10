from pathlib import Path

from crystal_transport.benchmark002 import run_benchmark


def test_small_benchmark002_writes_controls_and_figures(tmp_path: Path) -> None:
    result = run_benchmark(
        output_dir=tmp_path / "benchmark_002",
        resolution=6,
        crystalline_fractions=(0.0, 0.2),
        mobility_ratios=(1.0, 0.1),
        sweep_seeds=(7, 8),
        placement_seeds=(7, 8),
        grid_resolutions=(6,),
        contrast_matrix_values=(1.0e-3,),
        benchmark001_reference_dir=Path("results/benchmark_001"),
    )
    output = tmp_path / "benchmark_002"
    # The coarse smoke grid is an artifact-generation check.  The canonical
    # 32^3 B001 comparison is a release-gate control and is validated from the
    # committed full-resolution result, not from this platform-sensitive grid.
    assert result["validation"]["controls"]["same_seed_selection_digest_equal"]
    assert result["validation"]["controls"]["label_only_noop_field_exact"]
    assert (output / "manifest.json").exists()
    assert (output / "runs.csv").exists()
    assert (output / "summary.csv").exists()
    assert (output / "figures" / "mobility_sweep.png").exists()
    assert (output / "figures" / "placement_comparison.png").exists()
