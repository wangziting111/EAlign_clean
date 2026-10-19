#!/usr/bin/env python3
"""Paired-bootstrap comparison of methods in the pairwise evaluation.

The script reads saved per-pair scores; it does not rerun an alignment method.
With no arguments, it reproduces the paper comparisons of EAlign against
DeepAlign, TM-align, and DALI on HOMSTRAD, MALIDUP, and MALISAM.

Examples:

    python3 experiments/01_pairwise_alignment/analyze_pairwise_evaluation.py
    python3 experiments/01_pairwise_alignment/analyze_pairwise_evaluation.py \
        --method1 ealign --method2 deepalign --dataset HOMSTRAD

For every comparison, Delta F1 and its confidence interval are reported as
method1 minus method2 in percentage points. Empty F1 fields (tool errors) are
treated as zero, consistently with the all-pair accuracy protocol.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, Mapping, Sequence, Tuple

import numpy as np


REPO = Path(__file__).resolve().parents[2]
DEFAULT_DATASETS = ("homstrad", "malidup", "malisam")
DEFAULT_METHOD2 = ("deepalign", "tmalign", "dali")
METHOD_LABELS = {
    "ealign": "EAlign",
    "deepalign": "DeepAlign",
    "tmalign": "TM-align",
    "dali": "DALI",
}
REQUIRED_COLUMNS = {"benchmark", "method", "id", "f1"}


def method_name(value: str) -> str:
    normalized = value.strip().lower().replace("-", "").replace("_", "")
    aliases = {
        "ealign": "ealign",
        "deepalign": "deepalign",
        "tmalign": "tmalign",
        "dali": "dali",
        "dalilite": "dali",
    }
    if normalized not in aliases:
        raise argparse.ArgumentTypeError(
            f"unknown method {value!r}; choose EAlign, DeepAlign, TM-align, or DALI"
        )
    return aliases[normalized]


def dataset_name(value: str) -> str:
    normalized = value.strip().lower()
    if normalized not in DEFAULT_DATASETS:
        raise argparse.ArgumentTypeError(
            f"unknown dataset {value!r}; choose HOMSTRAD, MALIDUP, or MALISAM"
        )
    return normalized


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=REPO / "test_results/pairwise_evaluation/per_pair.csv",
    )
    parser.add_argument(
        "--method1",
        type=method_name,
        default="ealign",
        help="First method in method1 - method2 (default: EAlign).",
    )
    parser.add_argument(
        "--method2",
        action="append",
        type=method_name,
        default=[],
        help=(
            "Second method; repeat for multiple comparisons. The default is "
            "DeepAlign, TM-align, and DALI when method1 is EAlign."
        ),
    )
    parser.add_argument(
        "--dataset",
        action="append",
        type=dataset_name,
        default=[],
        help="Dataset to analyze; repeat as needed (default: all three).",
    )
    parser.add_argument("--bootstrap", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--confidence", type=float, default=0.95)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO / "test_results/pairwise_evaluation_statistics",
    )
    args = parser.parse_args()
    if args.bootstrap <= 0:
        parser.error("--bootstrap must be positive")
    if not 0.0 < args.confidence < 1.0:
        parser.error("--confidence must be between 0 and 1")
    if not args.method2:
        if args.method1 != "ealign":
            parser.error("--method2 is required when --method1 is not EAlign")
        args.method2 = list(DEFAULT_METHOD2)
    args.method2 = list(dict.fromkeys(args.method2))
    if args.method1 in args.method2:
        parser.error("--method1 and --method2 must differ")
    args.dataset = list(dict.fromkeys(args.dataset or DEFAULT_DATASETS))
    return args


def read_scores(path: Path) -> Dict[str, Dict[Tuple[str, str], float]]:
    if not path.is_file():
        raise SystemExit(f"per-pair CSV not found: {path}")
    scores: Dict[str, Dict[Tuple[str, str], float]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or ())
        if missing:
            raise SystemExit(f"{path} is missing columns: {sorted(missing)}")
        for line_number, row in enumerate(reader, start=2):
            method = row["method"].strip().lower()
            dataset = row["benchmark"].strip().lower()
            pair_id = row["id"].strip()
            key = (dataset, pair_id)
            method_scores = scores.setdefault(method, {})
            if key in method_scores:
                raise SystemExit(
                    f"duplicate row for {method}/{dataset}/{pair_id} "
                    f"in {path}:{line_number}"
                )
            try:
                method_scores[key] = float(row["f1"] or 0.0)
            except ValueError as exc:
                raise SystemExit(
                    f"invalid F1 for {method}/{dataset}/{pair_id} "
                    f"in {path}:{line_number}"
                ) from exc
    return scores


def paired_values(
    scores: Mapping[str, Mapping[Tuple[str, str], float]],
    method1: str,
    method2: str,
    dataset: str,
) -> Tuple[np.ndarray, np.ndarray]:
    for method in (method1, method2):
        if method not in scores:
            raise SystemExit(f"method {method!r} is absent from the input CSV")
    ids1 = sorted(pair_id for ds, pair_id in scores[method1] if ds == dataset)
    ids2 = sorted(pair_id for ds, pair_id in scores[method2] if ds == dataset)
    if not ids1:
        raise SystemExit(f"no rows for {method1}/{dataset}")
    if ids1 != ids2:
        missing = sorted(set(ids1) - set(ids2))
        extra = sorted(set(ids2) - set(ids1))
        raise SystemExit(
            f"pair mismatch for {dataset}, {method1} versus {method2}: "
            f"missing={missing[:5]} ({len(missing)} total), "
            f"extra={extra[:5]} ({len(extra)} total)"
        )
    values1 = np.array([scores[method1][(dataset, i)] for i in ids1])
    values2 = np.array([scores[method2][(dataset, i)] for i in ids1])
    return values1, values2


def paired_bootstrap_ci(
    differences: np.ndarray,
    resamples: int,
    confidence: float,
    seed: int,
) -> Tuple[float, float]:
    """Percentile CI of the mean paired difference, in percentage points."""
    rng = np.random.default_rng(seed)
    means = np.empty(resamples, dtype=float)
    batch_size = min(1_000, resamples)
    for start in range(0, resamples, batch_size):
        count = min(batch_size, resamples - start)
        indices = rng.integers(
            0, len(differences), size=(count, len(differences))
        )
        means[start : start + count] = differences[indices].mean(axis=1)
    tail = (1.0 - confidence) / 2.0
    low, high = np.quantile(means, (tail, 1.0 - tail))
    return float(low), float(high)


def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    scores = read_scores(args.input_csv)
    rows = []
    for dataset in args.dataset:
        for method2 in args.method2:
            values1, values2 = paired_values(
                scores, args.method1, method2, dataset
            )
            differences = 100.0 * (values1 - values2)
            low, high = paired_bootstrap_ci(
                differences, args.bootstrap, args.confidence, args.seed
            )
            rows.append(
                {
                    "dataset": dataset,
                    "method1": args.method1,
                    "method2": method2,
                    "comparison": (
                        f"{METHOD_LABELS[args.method1]} - {METHOD_LABELS[method2]}"
                    ),
                    "pairs": len(differences),
                    "method1_macro_f1_pct": 100.0 * float(values1.mean()),
                    "method2_macro_f1_pct": 100.0 * float(values2.mean()),
                    "delta_f1_pct_points": float(differences.mean()),
                    "ci_low_pct_points": low,
                    "ci_high_pct_points": high,
                }
            )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "pairwise_paired_bootstrap_ci.csv"
    json_path = args.output_dir / "pairwise_paired_bootstrap_ci.json"
    write_csv(csv_path, rows)
    json_path.write_text(
        json.dumps(
            {
                "definition": (
                    "Delta F1 = macro per-pair F1(method1) - "
                    "macro per-pair F1(method2), in percentage points"
                ),
                "ci": "percentile paired-bootstrap CI of the mean paired difference",
                "empty_f1": "treated as zero, including tool errors",
                "bootstrap": {
                    "resamples": args.bootstrap,
                    "seed": args.seed,
                    "confidence": args.confidence,
                    "rng": "numpy.random.default_rng (PCG64)",
                },
                "input": str(args.input_csv.resolve()),
                "comparisons": rows,
            },
            indent=2,
        )
        + "\n"
    )

    print("Dataset   Comparison                N   Delta F1 (pp)   95% CI (pp)")
    for row in rows:
        print(
            f"{row['dataset'].upper():<9} "
            f"{row['comparison']:<25} "
            f"{row['pairs']:4d}   "
            f"{row['delta_f1_pct_points']:+12.3f}   "
            f"[{row['ci_low_pct_points']:.3f}, "
            f"{row['ci_high_pct_points']:.3f}]"
        )
    print(f"\nWrote {csv_path}")
    print(f"Wrote {json_path}")


if __name__ == "__main__":
    main()
