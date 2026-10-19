#!/usr/bin/env python3
"""Score saved residue correspondences with the native EAlign objective.

For every structure pair, the script reads the saved EAlign, DeepAlign,
TM-align, and DALI correspondences.  Each method's correspondences are refit
with the same least-squares Kabsch transform, then scored by

    E_h(i,j|T) = W(i,j) V(i,j|R) g(d_ij(T)).

The C++ helper reuses EAlign's own 3Di, substitution, orientation, distance,
and pairwise min-max normalization code.  No alignment method is rerun.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from statistics import mean

from evaluate_pairwise_alignment import binary_prediction, dali_prediction
from run_pairwise_alignment import BENCHMARKS, EXPECTED, REPO, Case, load_cases


METHODS = ("ealign", "deepalign", "tmalign", "dali")
METHOD_LABELS = {
    "ealign": "EAlign",
    "deepalign": "DeepAlign",
    "tmalign": "TM-align",
    "dali": "DALI",
}
BINARY_ROOTS = {
    method: REPO / "test_results" / "pairwise_runs"
    for method in ("ealign", "deepalign", "tmalign")
}
DALI_ROOTS = {
    "malisam": REPO / "test_results" / "malisam_dali",
    "malidup": REPO / "test_results" / "malidup_dali",
    "homstrad": REPO / "test_results" / "homstrad_dali",
}


def load_prediction(case: Case, method: str) -> tuple[str, set[tuple[int, int]]]:
    if method == "dali":
        return dali_prediction(case, DALI_ROOTS[case.benchmark])
    return binary_prediction(case, method, BINARY_ROOTS[method])


def parse_scores(text: str) -> dict[str, dict[str, str]]:
    lines = text.splitlines()
    header = next(
        (index for index, line in enumerate(lines) if line.startswith("method\tlen1\t")),
        None,
    )
    if header is None:
        raise RuntimeError("scorer output has no score-table header")
    rows = list(csv.DictReader(lines[header:], delimiter="\t"))
    return {row["method"]: row for row in rows if row.get("len1")}


def score_case(
    case: Case,
    methods: tuple[str, ...],
    scorer: Path,
    timeout: int,
) -> list[dict]:
    base = {
        "dataset": case.benchmark,
        "case_id": case.case_id,
        "pdb_a": str(case.pdb_a),
        "pdb_b": str(case.pdb_b),
    }
    predictions: dict[str, set[tuple[int, int]]] = {}
    rows: dict[str, dict] = {}
    for method in methods:
        try:
            status, pairs = load_prediction(case, method)
            rows[method] = {**base, "method": method, "status": status, "error": ""}
            if status in ("ok", "hit"):
                if not pairs:
                    raise ValueError("successful alignment contains no residue pairs")
                predictions[method] = pairs
        except (OSError, ValueError, KeyError, IndexError) as error:
            rows[method] = {
                **base,
                "method": method,
                "status": "evaluation_error",
                "error": str(error),
            }

    if predictions:
        with tempfile.TemporaryDirectory(prefix=f"ealign_score_{case.case_id}_") as tmp:
            pair_file = Path(tmp) / "pairs.tsv"
            with pair_file.open("w") as handle:
                handle.write("method\tquery_index\ttarget_index\n")
                for method in methods:
                    for query, target in sorted(predictions.get(method, set())):
                        handle.write(f"{method}\t{query}\t{target}\n")
            command = [
                str(scorer),
                "--score-pairs",
                str(case.pdb_a),
                str(case.pdb_b),
                str(pair_file),
            ]
            try:
                done = subprocess.run(
                    command, capture_output=True, text=True, timeout=timeout
                )
                if done.returncode != 0:
                    raise RuntimeError(
                        f"scorer exit {done.returncode}: {done.stderr.strip()}"
                    )
                scored = parse_scores(done.stdout)
                for method, pairs in predictions.items():
                    score = scored.get(method)
                    if score is None:
                        raise RuntimeError(f"scorer omitted {method}")
                    if int(score["aligned_pairs"]) != len(pairs):
                        raise RuntimeError(
                            f"{method} pair-count mismatch: "
                            f"{score['aligned_pairs']} != {len(pairs)}"
                        )
                    rows[method].update(score)
            except (OSError, RuntimeError, subprocess.TimeoutExpired) as error:
                for method in predictions:
                    rows[method].update(status="scoring_error", error=str(error))

    columns = (
        "len1",
        "len2",
        "aligned_pairs",
        "score_sum",
        "mean_per_aligned_pair",
        "score_per_shorter_length",
    )
    return [
        {**rows[method], **{key: rows[method].get(key, "") for key in columns}}
        for method in methods
    ]


def summarize(rows: list[dict], benchmarks: tuple[str, ...], methods: tuple[str, ...]) -> list[dict]:
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["dataset"], row["method"])].append(row)

    summaries = []
    for dataset in benchmarks:
        for method in methods:
            selected = grouped[(dataset, method)]
            measured = [row for row in selected if row["mean_per_aligned_pair"] != ""]
            if len(selected) != EXPECTED[dataset]:
                raise ValueError(
                    f"Expected {EXPECTED[dataset]} {dataset}/{method} rows, "
                    f"got {len(selected)}"
                )
            statuses = Counter(row["status"] for row in selected)
            summaries.append(
                {
                    "dataset": dataset,
                    "method": method,
                    "method_label": METHOD_LABELS[method],
                    "requested_pairs": len(selected),
                    "aligned_pairs": len(measured),
                    "coverage": len(measured) / len(selected),
                    "ok": statuses["ok"],
                    "hit": statuses["hit"],
                    "no_hit": statuses["no_hit"],
                    "error": statuses["error"],
                    "timeout": statuses["timeout"],
                    "evaluation_error": statuses["evaluation_error"],
                    "scoring_error": statuses["scoring_error"],
                    "mean_per_aligned_pair_score_hit_only": mean(
                        float(row["mean_per_aligned_pair"]) for row in measured
                    ),
                    "pooled_per_residue_pair_score_hit_only": sum(
                        float(row["score_sum"]) for row in measured
                    )
                    / sum(int(row["aligned_pairs"]) for row in measured),
                }
            )
    return summaries


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def print_summary(rows: list[dict]) -> None:
    print("Dataset   Method      scored/N   Mean per-aligned-pair EAlign score")
    for row in rows:
        print(
            f"{row['dataset']:<9} {row['method_label']:<11} "
            f"{row['aligned_pairs']:>4}/{row['requested_pairs']:<4} "
            f"{row['mean_per_aligned_pair_score_hit_only']:>34.4f}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--benchmarks", nargs="+", choices=BENCHMARKS, default=list(BENCHMARKS)
    )
    parser.add_argument("--methods", nargs="+", choices=METHODS, default=list(METHODS))
    parser.add_argument(
        "--scorer", type=Path, default=REPO / "binary" / "ealign_score_pairs"
    )
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO / "test_results" / "ealign_alignment_score",
    )
    args = parser.parse_args()
    scorer = args.scorer.resolve()
    if not scorer.is_file():
        raise FileNotFoundError(f"Build the scorer first: make score-pairs ({scorer})")

    benchmarks = tuple(args.benchmarks)
    methods = tuple(args.methods)
    cases = load_cases(list(benchmarks))
    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(score_case, case, methods, scorer, args.timeout): case
            for case in cases
        }
        for completed, future in enumerate(as_completed(futures), 1):
            rows.extend(future.result())
            if completed % 100 == 0 or completed == len(futures):
                print(f"scored {completed}/{len(futures)} structure pairs", flush=True)

    rows.sort(key=lambda row: (row["dataset"], row["method"], row["case_id"]))
    summaries = summarize(rows, benchmarks, methods)
    output = args.output_dir.resolve()
    write_csv(output / "per_pair.csv", rows)
    write_csv(output / "summary.csv", summaries)
    print_summary(summaries)
    print(f"per-pair: {output / 'per_pair.csv'}")
    print(f"summary:  {output / 'summary.csv'}")


if __name__ == "__main__":
    main()
