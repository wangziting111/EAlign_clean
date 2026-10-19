#!/usr/bin/env python3
"""Build and optionally benchmark EAlign score ablations.

Examples:
  python3 experiments/02_score_ablation/ablation_ealign_score.py --all --build-only
  python3 experiments/02_score_ablation/ablation_ealign_score.py --without-aa --workers 8
  python3 experiments/02_score_ablation/ablation_ealign_score.py --without-vector --limit 20

Without a setting flag, all four ablations are selected. Unless --build-only is
used, each selected binary is evaluated by run_alignment_benchmark.py.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
RUNNER = Path(__file__).resolve().with_name("run_alignment_benchmark.py")


@dataclass(frozen=True)
class Setting:
    name: str
    binary_name: str
    evidence_mode: int
    geometry_mode: int


SETTINGS = {
    "baseline": Setting("baseline", "ealign", 0, -1),
    "without-aa": Setting("without-aa", "ealign_without_aa", 1, -1),
    "without-3di": Setting("without-3di", "ealign_without_3di", 2, -1),
     "without-3di-aa": Setting("without-3di-aa", "ealign_without_3di_aa", 4, -1),
    "without-vector": Setting(
        "without-vector", "ealign_without_vector", 0, 0
    ),
    "without-distance": Setting(
        "without-distance", "ealign_without_distance", 0, 1
    ),
}
ABLATIONS = tuple(name for name in SETTINGS if name != "baseline")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--setting",
        action="append",
        choices=tuple(SETTINGS),
        help="Select a setting; repeat this option to select several.",
    )
    parser.add_argument("--baseline", action="store_true")
    parser.add_argument("--without-aa", action="store_true")
    parser.add_argument("--without-3di", action="store_true")
    parser.add_argument("--without-3di-aa", action="store_true")
    parser.add_argument("--without-vector", action="store_true")
    parser.add_argument(
        "--without-distance",
        "--without-g",
        dest="without_distance",
        action="store_true",
    )
    parser.add_argument("--all", action="store_true", help="Select all four ablations.")
    parser.add_argument(
        "--build-only",
        action="store_true",
        help="Build binaries without running the full alignment benchmark.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=max(1, min(8, os.cpu_count() or 1)),
        help="Parallel compiler jobs (default: up to 8).",
    )
    parser.add_argument("--workers", type=int, default=8, help="Benchmark workers.")
    parser.add_argument(
        "--benchmarks",
        nargs="+",
        choices=("malisam", "malidup", "homstrad"),
        default=("malisam", "malidup", "homstrad"),
    )
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument(
        "--reuse-existing-results",
        action="store_true",
        help="Allow the benchmark runner to reuse compatible cached outputs.",
    )
    return parser.parse_args()


def selected_settings(args: argparse.Namespace) -> list[Setting]:
    names = list(args.setting or [])
    for enabled, name in (
        (args.baseline, "baseline"),
        (args.without_aa, "without-aa"),
        (args.without_3di, "without-3di"),
        (args.without_3di_aa, "without-3di-aa"),
        (args.without_vector, "without-vector"),
        (args.without_distance, "without-distance"),
    ):
        if enabled:
            names.append(name)
    if args.all or not names:
        names.extend(ABLATIONS)
    return [SETTINGS[name] for name in dict.fromkeys(names)]


def run(command: list[str]) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=REPO, check=True)


def build(setting: Setting, jobs: int) -> Path:
    binary_dir = REPO / "binary"
    binary_dir.mkdir(exist_ok=True)
    final_binary = binary_dir / setting.binary_name
    temporary_target = f"binary/.{setting.binary_name}.building-{os.getpid()}"

    run(["make", "clean-objects"])
    run(
        [
            "make",
            f"-j{jobs}",
            f"TARGET={temporary_target}",
            f"EALIGN_ALIGNMENT_EVIDENCE_MODE={setting.evidence_mode}",
            f"EALIGN_ALIGNMENT_GEOMETRY_MODE={setting.geometry_mode}",
        ]
    )
    (REPO / temporary_target).replace(final_binary)
    print(f"built {setting.name}: {final_binary}", flush=True)
    return final_binary


def benchmark(binary: Path, args: argparse.Namespace) -> None:
    command = [
        sys.executable,
        str(RUNNER),
        "--binary",
        str(binary),
        "--output-root",
        str(REPO / "test_results" / f"full_alignment_{binary.name}"),
        "--workers",
        str(args.workers),
        "--benchmarks",
        *args.benchmarks,
    ]
    if args.limit:
        command.extend(("--limit", str(args.limit)))
    if not args.reuse_existing_results:
        command.append("--overwrite")
    run(command)


def main() -> None:
    args = parse_args()
    if args.jobs < 1 or args.workers < 1 or args.limit < 0:
        raise SystemExit("--jobs/--workers must be positive and --limit non-negative")
    if not RUNNER.is_file():
        raise SystemExit(f"benchmark runner not found: {RUNNER}")

    for setting in selected_settings(args):
        binary = build(setting, args.jobs)
        if not args.build_only:
            benchmark(binary, args)


if __name__ == "__main__":
    main()
