#!/usr/bin/env python3
"""Rescore saved pairwise alignments with the normalized EAlign objective.

The script reads residue correspondences from ``test_results/pairwise_runs``,
refits each correspondence set with Kabsch, and evaluates

    EAlignScore = sum(normalized(BLOSUM + 3Di) * orientation * distance).

``per_pair.csv`` reports the score sum, its mean over aligned residue pairs,
and the sum normalized by the shorter structure length. ``summary.csv``
reports dataset-level macro means. Failed/no-hit alignments remain in the
per-pair output but are excluded from score means.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from statistics import mean

from evaluate_pairwise_alignment import (
    DEFAULT_METHODS,
    binary_prediction,
    dali_prediction,
)
from run_pairwise_alignment import BENCHMARKS, EXPECTED, REPO, Case, load_cases


METHODS = DEFAULT_METHODS
METHOD_LABELS = {
    "ealign": "EAlign",
    "deepalign": "DeepAlign",
    "tmalign": "TM-align",
    "dali": "DALI",
    "gtalign_cuda": "GTalign-CUDA",
    "gtalign_patched_cpu": "GTalign-CPU-1",
    "gtalign_patched_cpu_t20": "GTalign-CPU-20",
    "foldseek": "Foldseek",
}
SUCCESS = {"ok", "hit"}
PER_PAIR_FIELDS = (
    "dataset",
    "case_id",
    "method",
    "method_label",
    "status",
    "len1",
    "len2",
    "aligned_pairs",
    "score_sum",
    "mean_per_aligned_pair",
    "mean_normalized_shorter",
    "error",
)


def load_prediction(
    case: Case, method: str, pairwise_root: Path
) -> tuple[str, set[tuple[int, int]]]:
    if method == "dali":
        return dali_prediction(case, pairwise_root / "dali" / case.benchmark)
    return binary_prediction(case, method, pairwise_root)


def parse_score_table(stdout: str) -> dict[str, dict[str, str]]:
    lines = stdout.splitlines()
    start = next(
        (i for i, line in enumerate(lines) if line.startswith("method\tlen1\t")),
        None,
    )
    if start is None:
        raise RuntimeError("EAlign scorer returned no score table")
    return {
        row["method"]: row
        for row in csv.DictReader(lines[start:], delimiter="\t")
        if row.get("method")
    }


def run_scorer(
    case: Case,
    predictions: dict[str, set[tuple[int, int]]],
    scorer: Path,
    timeout: int,
) -> dict[str, dict[str, str]]:
    with tempfile.TemporaryDirectory(prefix=f"ealign_score_{case.case_id}_") as tmp:
        pair_file = Path(tmp) / "pairs.tsv"
        with pair_file.open("w", encoding="utf-8") as handle:
            handle.write("method\tquery_index\ttarget_index\n")
            for method, pairs in predictions.items():
                for query, target in sorted(pairs):
                    handle.write(f"{method}\t{query}\t{target}\n")

        completed = subprocess.run(
            [
                str(scorer),
                "--score-pairs",
                str(case.pdb_a),
                str(case.pdb_b),
                str(pair_file),
            ],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    if completed.returncode != 0:
        raise RuntimeError(
            f"scorer exit {completed.returncode}: {completed.stderr.strip()}"
        )
    return parse_score_table(completed.stdout)


def score_case(
    case: Case,
    methods: tuple[str, ...],
    pairwise_root: Path,
    scorer: Path,
    timeout: int,
) -> list[dict]:
    rows = {}
    predictions = {}
    for method in methods:
        row = {
            "dataset": case.benchmark,
            "case_id": case.case_id,
            "method": method,
            "method_label": METHOD_LABELS[method],
            "status": "",
            "len1": "",
            "len2": "",
            "aligned_pairs": "",
            "score_sum": "",
            "mean_per_aligned_pair": "",
            "mean_normalized_shorter": "",
            "error": "",
        }
        try:
            status, pairs = load_prediction(case, method, pairwise_root)
            row["status"] = status
            if status in SUCCESS:
                if len(pairs) < 3:
                    raise ValueError("successful alignment contains fewer than three pairs")
                predictions[method] = pairs
        except (OSError, ValueError, KeyError, IndexError) as error:
            row.update(status="evaluation_error", error=str(error))
        rows[method] = row

    if predictions:
        try:
            scored = run_scorer(case, predictions, scorer, timeout)
            for method, pairs in predictions.items():
                values = scored.get(method)
                if values is None:
                    raise RuntimeError(f"scorer omitted {method}")
                aligned_pairs = int(values["aligned_pairs"])
                if aligned_pairs != len(pairs):
                    raise RuntimeError(
                        f"{method} pair-count mismatch: {aligned_pairs} != {len(pairs)}"
                    )
                len1 = int(values["len1"])
                len2 = int(values["len2"])
                score_sum = float(values["score_sum"])
                rows[method].update(
                    len1=len1,
                    len2=len2,
                    aligned_pairs=aligned_pairs,
                    score_sum=score_sum,
                    mean_per_aligned_pair=float(values["mean_per_aligned_pair"]),
                    mean_normalized_shorter=score_sum / min(len1, len2),
                )
        except (OSError, RuntimeError, subprocess.TimeoutExpired) as error:
            for method in predictions:
                rows[method].update(status="scoring_error", error=str(error))

    return [rows[method] for method in methods]


def summarize(
    rows: list[dict], benchmarks: tuple[str, ...], methods: tuple[str, ...]
) -> list[dict]:
    summary = []
    for dataset in benchmarks:
        for method in methods:
            selected = [
                row
                for row in rows
                if row["dataset"] == dataset and row["method"] == method
            ]
            if len(selected) != EXPECTED[dataset]:
                raise ValueError(
                    f"Expected {EXPECTED[dataset]} {dataset}/{method} rows, "
                    f"found {len(selected)}"
                )
            scored = [row for row in selected if row["score_sum"] != ""]
            summary.append(
                {
                    "dataset": dataset,
                    "method": method,
                    "method_label": METHOD_LABELS[method],
                    "total_cases": len(selected),
                    "scored_cases": len(scored),
                    "mean_score_sum": mean(row["score_sum"] for row in scored),
                    # This macro mean is the EAlign-score column in Table 1.
                    "mean_per_aligned_pair": mean(
                        row["mean_per_aligned_pair"] for row in scored
                    ),
                    "mean_normalized_shorter": mean(
                        row["mean_normalized_shorter"] for row in scored
                    ),
                }
            )
    return summary


def write_csv(path: Path, rows: list[dict], fields: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def print_summary(rows: list[dict]) -> None:
    print(
        "Dataset   Method          scored/total  mean(score_sum)  "
        "mean/pair  mean/shorter"
    )
    for row in rows:
        print(
            f"{row['dataset']:<9} {row['method_label']:<15} "
            f"{row['scored_cases']:>4}/{row['total_cases']:<4} "
            f"{row['mean_score_sum']:>16.6f} "
            f"{row['mean_per_aligned_pair']:>10.6f} "
            f"{row['mean_normalized_shorter']:>13.6f}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--benchmarks", nargs="+", choices=BENCHMARKS, default=list(BENCHMARKS)
    )
    parser.add_argument(
        "--methods", nargs="+", choices=METHODS, default=list(METHODS)
    )
    parser.add_argument(
        "--pairwise-root",
        type=Path,
        default=REPO / "test_results/pairwise_runs",
    )
    parser.add_argument("--scorer", type=Path, default=REPO / "binary/ealign")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO / "test_results/01_pairwise_alignment/ealign_score",
    )
    args = parser.parse_args()

    pairwise_root = args.pairwise_root.resolve()
    scorer = args.scorer.resolve()
    if not pairwise_root.is_dir():
        raise FileNotFoundError(pairwise_root)
    if not scorer.is_file():
        raise FileNotFoundError(f"Build EAlign first with make: {scorer}")

    benchmarks = tuple(args.benchmarks)
    methods = tuple(args.methods)
    cases = load_cases(list(benchmarks))
    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [
            pool.submit(
                score_case,
                case,
                methods,
                pairwise_root,
                scorer,
                args.timeout,
            )
            for case in cases
        ]
        for completed, future in enumerate(as_completed(futures), 1):
            rows.extend(future.result())
            if completed % 100 == 0 or completed == len(futures):
                print(f"scored {completed}/{len(futures)} structure pairs", flush=True)

    rows.sort(key=lambda row: (row["dataset"], row["method"], row["case_id"]))
    summary = summarize(rows, benchmarks, methods)
    output = args.output_dir.resolve()
    write_csv(output / "per_pair.csv", rows, PER_PAIR_FIELDS)
    write_csv(output / "summary.csv", summary, tuple(summary[0]))
    print_summary(summary)
    print(f"per-pair: {output / 'per_pair.csv'}")
    print(f"summary:  {output / 'summary.csv'}")


if __name__ == "__main__":
    main()
