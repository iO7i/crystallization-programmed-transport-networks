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
    assert result["validation"]["controls"]["benchmark001_baseline_pass"]
    assert result["validation"]["controls"]["d_crystal_equals_d_mobile_noop_pass"]
    assert (output / "manifest.json").exists()
    assert (output / "runs.csv").exists()
    assert (output / "summary.csv").exists()
    assert (output / "figures" / "mobility_sweep.png").exists()
    assert (output / "figures" / "placement_comparison.png").exists()
