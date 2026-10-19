#!/usr/bin/env python3
"""Summarize single-worker pairwise runtime experiments by benchmark."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from statistics import mean, median

from run_pairwise_alignment import BENCHMARKS, REPO, TIMING_SCOPE, ca_sequence, load_cases

METHODS = ("ealign", "deepalign", "tmalign", "dali")
DEFAULT_ROOTS = {
    "ealign": REPO / "test_results/runtime_ealign",
    "deepalign": REPO / "test_results/runtime_deepalign",
    "tmalign": REPO / "test_results/runtime_tmalign",
}
DEFAULT_DALI_ROOTS = {
    benchmark: REPO / (
        f"test_results/{benchmark}_dali"
    )
    for benchmark in BENCHMARKS
}
DEFAULT_OUTPUT = REPO / "test_results/runtime_summary.csv"


def dataset_statistics(benchmarks: list[str], split_root: Path) -> dict[str, dict]:
    cases = load_cases(benchmarks, split_root)
    result = {}
    for benchmark in benchmarks:
        selected = [case for case in cases if case.benchmark == benchmark]
        pair_lengths = [
            (len(ca_sequence(case.pdb_a)) + len(ca_sequence(case.pdb_b))) / 2
            for case in selected
        ]
        result[benchmark] = {
            "cases": {case.case_id: case for case in selected},
            "dataset_size": len(selected),
            "mean_pair_length": mean(pair_lengths),
        }
    return result


def load_runtime_rows(root: Path, method: str, benchmark: str) -> dict[str, dict]:
    rows = {}
    pattern = "pairs/*/result.json" if method == "dali" else f"{method}/{benchmark}/*/result.json"
    for path in sorted(root.glob(pattern)):
        result = json.loads(path.read_text())
        saved_method = "dali" if method == "dali" else result.get("method")
        saved_benchmark = result.get("benchmark", "homstrad" if method == "dali" else None)
        if saved_method != method or saved_benchmark != benchmark:
            raise ValueError(f"Metadata does not match its directory: {path}")
        case_id = result.get("family") if method == "dali" else result.get("id")
        if not case_id or case_id in rows:
            raise ValueError(f"Missing or duplicate case ID in {path}")
        rows[case_id] = result
    return rows


def validated_runtimes(
    method: str,
    root: Path,
    benchmark: str,
    dataset: dict,
    allow_incomplete: bool,
) -> tuple[dict[str, float], list[str], list[str]]:
    saved = load_runtime_rows(root, method, benchmark)
    expected_ids = set(dataset["cases"])
    extra = sorted(set(saved) - expected_ids)
    if extra:
        raise ValueError(f"Unexpected {benchmark}/{method} cases: {extra[:10]}")

    successful = {}
    failures = []
    for case_id, result in saved.items():
        if method != "dali" and result.get("timing_scope") != TIMING_SCOPE:
            raise ValueError(
                f"Wrong timing scope for {benchmark}/{method}/{case_id}: "
                f"{result.get('timing_scope')!r}"
            )
        case = dataset["cases"][case_id]
        expected_paths = (case.pdb_a.resolve(), case.pdb_b.resolve())
        saved_paths = (Path(result["pdb_a"]).resolve(), Path(result["pdb_b"]).resolve())
        if saved_paths != expected_paths and not (method == "dali" and saved_paths == expected_paths[::-1]):
            raise ValueError(f"Input pair changed for {benchmark}/{method}/{case_id}")
        valid_statuses = {"hit", "no_hit"} if method == "dali" else {"ok"}
        if result.get("status") in valid_statuses:
            seconds = float(result["seconds"])
            if seconds < 0:
                raise ValueError(f"Negative runtime for {benchmark}/{method}/{case_id}")
            successful[case_id] = seconds
        else:
            failures.append(case_id)

    missing = sorted(expected_ids - set(saved))
    # Dali import/execution failures are a benchmark result: report and exclude
    # them from runtime averages. Missing result files still mean an incomplete run.
    blocking_failures = failures if method != "dali" else []
    if (missing or blocking_failures) and not allow_incomplete:
        raise RuntimeError(
            f"{benchmark}/{method} is incomplete: "
            f"missing={len(missing)}, failed={len(failures)}. "
            "Wait for the run to finish, or use --allow-incomplete for a preview."
        )
    return successful, failures, missing


def summarize_method(
    method: str,
    root: Path,
    benchmark: str,
    dataset: dict,
    allow_incomplete: bool,
) -> dict:
    successful, failures, missing = validated_runtimes(
        method, root, benchmark, dataset, allow_incomplete
    )
    values = list(successful.values())
    return {
        "benchmark": benchmark,
        "method": method,
        "dataset_size": dataset["dataset_size"],
        "measured_pairs": len(values),
        "failed_pairs": len(failures),
        "missing_pairs": len(missing),
        "mean_pair_length": dataset["mean_pair_length"],
        "mean_runtime_seconds": mean(values) if values else None,
        "median_runtime_seconds": median(values) if values else None,
        "timing_scope": (
            "dali.pl_subprocess_wall_time" if method == "dali" else TIMING_SCOPE
        ),
    }


def show(rows: list[dict], benchmarks: list[str]) -> None:
    for benchmark in benchmarks:
        selected = [row for row in rows if row["benchmark"] == benchmark]
        first = selected[0]
        print(f"\n[{benchmark}]")
        print(f"dataset size: {first['dataset_size']}")
        print(f"mean (L1 + L2) / 2: {first['mean_pair_length']:.3f}")
        print("method       measured  failed  mean runtime (s)  median runtime (s)")
        for row in selected:
            mean_value = row["mean_runtime_seconds"]
            median_value = row["median_runtime_seconds"]
            print(
                f"{row['method']:<12} {row['measured_pairs']:>8}  {row['failed_pairs']:>6}  "
                f"{mean_value if mean_value is not None else float('nan'):>16.6f}  "
                f"{median_value if median_value is not None else float('nan'):>18.6f}"
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmarks", nargs="+", choices=BENCHMARKS, default=list(BENCHMARKS))
    parser.add_argument("--ealign-root", type=Path, default=DEFAULT_ROOTS["ealign"])
    parser.add_argument("--deepalign-root", type=Path, default=DEFAULT_ROOTS["deepalign"])
    parser.add_argument("--tmalign-root", type=Path, default=DEFAULT_ROOTS["tmalign"])
    parser.add_argument("--dali-malisam-root", type=Path, default=DEFAULT_DALI_ROOTS["malisam"])
    parser.add_argument("--dali-malidup-root", type=Path, default=DEFAULT_DALI_ROOTS["malidup"])
    parser.add_argument("--dali-homstrad-root", type=Path, default=DEFAULT_DALI_ROOTS["homstrad"])
    parser.add_argument(
        "--split-root",
        type=Path,
        default=REPO / "test_results/HOMSTRAD_split_structures",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help="Preview currently completed pairs instead of requiring a complete run",
    )
    args = parser.parse_args()

    datasets = dataset_statistics(args.benchmarks, args.split_root)
    roots = {
        method: getattr(args, f"{method}_root").resolve()
        for method in METHODS if method != "dali"
    }
    dali_roots = {
        benchmark: getattr(args, f"dali_{benchmark}_root").resolve()
        for benchmark in BENCHMARKS
    }
    rows = [
        summarize_method(
                method,
                dali_roots[benchmark] if method == "dali" else roots[method],
                benchmark,
                datasets[benchmark],
                args.allow_incomplete,
            )
        for benchmark in args.benchmarks
        for method in METHODS
    ]

    args.output = args.output.resolve()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = (
        "benchmark", "method", "dataset_size", "measured_pairs", "failed_pairs",
        "missing_pairs", "mean_pair_length", "mean_runtime_seconds",
        "median_runtime_seconds", "timing_scope",
    )
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    show(rows, args.benchmarks)
    print(f"\nCSV: {args.output}")


if __name__ == "__main__":
    main()
