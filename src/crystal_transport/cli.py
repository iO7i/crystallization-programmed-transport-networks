"""Command-line entry point for the reproducible benchmark."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .benchmark001 import run_benchmark
from .benchmark002 import run_benchmark as run_benchmark002
from .validation import run_validation, validation_passes


def main() -> int:
    parser = argparse.ArgumentParser(prog="crystal-transport")
    subparsers = parser.add_subparsers(dest="command", required=True)
    validate_parser = subparsers.add_parser("validate", help="run analytical solver validations")
    validate_parser.add_argument("--resolution", type=int, default=16)
    benchmark_parser = subparsers.add_parser("benchmark", help="run a named benchmark")
    benchmark_subparsers = benchmark_parser.add_subparsers(dest="benchmark_command", required=True)
    run_parser = benchmark_subparsers.add_parser("run", help="execute a benchmark")
    run_parser.add_argument("name", choices=["benchmark_001", "benchmark_002"])
    run_parser.add_argument("--output", type=Path, default=None)
    run_parser.add_argument("--resolution", type=int, default=32)
    run_parser.add_argument(
        "--validation-profile",
        choices=["canonical", "smoke"],
        default="canonical",
        help="Benchmark 002 acceptance profile; smoke is for fast CI artifact checks",
    )
    args = parser.parse_args()
    if args.command == "validate":
        result = run_validation(shape=(args.resolution,) * 3)
        print(json.dumps(result, indent=2))
        return 0 if validation_passes(result) else 1
    output = args.output or Path("results") / args.name
    if args.name == "benchmark_001":
        result = run_benchmark(output_dir=output, resolution=args.resolution)
    else:
        result = run_benchmark002(
            output_dir=output,
            resolution=args.resolution,
            validation_profile=args.validation_profile,
        )
    print(
        json.dumps(
            {"validation_passed": result["validation_passed"], "output": str(args.output)},
            indent=2,
        )
    )
    return 0 if result["validation_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
