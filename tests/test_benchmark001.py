from crystal_transport.benchmark001 import run_benchmark


def test_small_benchmark_writes_reproducible_artifacts(tmp_path) -> None:
    result = run_benchmark(
        output_dir=tmp_path / "benchmark_001",
        resolution=8,
        convergence_resolutions=(8,),
    )
    assert result["validation_passed"]
    assert (tmp_path / "benchmark_001" / "manifest.json").exists()
    assert (tmp_path / "benchmark_001" / "figures" / "benchmark_001_overview.png").exists()

