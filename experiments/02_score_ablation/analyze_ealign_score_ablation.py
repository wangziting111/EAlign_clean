#!/usr/bin/env python3
"""Reproduce the EAlign-score ablation table and paired-bootstrap CIs.

Run from the EAlign repository root:

    python3 experiments/02_score_ablation/analyze_ealign_score_ablation.py

The script reads the saved per-pair benchmark CSVs; it does not rerun EAlign.
For every dataset and ablation, it requires exactly the same successful pair
IDs as Full EAlign.  The evaluation reference is restricted to alignment
columns in which both released sequences are uppercase.  It reports uppercase
macro precision, recall, and F1, the signed change
Delta F1 = F1(variant) - F1(full), the number of pairs with uppercase F1 < 0.5,
and a percentile paired-bootstrap 95% CI for
F1(full) - F1(variant).  All F1 differences are percentage points.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence, Tuple

import numpy as np


REPO = Path(__file__).resolve().parents[2]
DEFAULT_VARIANTS = (
    (
        "without_aa",
        REPO
        / "test_results/full_alignment_ealign_without_aa"
        / "ealign_without_aa_alignment_benchmark.csv",
    ),
    (
        "without_3di",
        REPO
        / "test_results/full_alignment_ealign_without_3di"
        / "ealign_without_3di_alignment_benchmark.csv",
    ),
    (
        "without_vector",
        REPO
        / "test_results/full_alignment_ealign_without_vector"
        / "ealign_without_vector_alignment_benchmark.csv",
    ),
    (
        "without_distance",
        REPO
        / "test_results/full_alignment_ealign_without_distance"
        / "ealign_without_distance_alignment_benchmark.csv",
    ),
)
DEFAULT_DATASETS = ("homstrad", "malidup", "malisam")
METRIC_COLUMNS = {
    "precision": "uppercase_precision",
    "recall": "uppercase_recall",
    "f1": "uppercase_f1",
}
REQUIRED_COLUMNS = {"benchmark", "id", "status", *METRIC_COLUMNS.values()}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full-csv",
        type=Path,
        default=REPO
        / "test_results/full_alignment_ealign/ealign_alignment_benchmark.csv",
    )
    parser.add_argument(
        "--variant",
        action="append",
        default=[],
        metavar="LABEL=CSV",
        help=(
            "Variant label and benchmark CSV. Repeat to replace the four "
            "paper defaults."
        ),
    )
    parser.add_argument(
        "--dataset",
        action="append",
        default=[],
        help="Dataset to include; repeat as needed (default: paper datasets).",
    )
    parser.add_argument("--bootstrap", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--confidence", type=float, default=0.95)
    parser.add_argument("--bad-threshold", type=float, default=0.5)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO / "test_results/ealign_score_ablation_statistics",
    )
    args = parser.parse_args()
    if args.bootstrap <= 0:
        parser.error("--bootstrap must be positive")
    if not 0.0 < args.confidence < 1.0:
        parser.error("--confidence must be between 0 and 1")
    return args


def parse_variants(values: Sequence[str]) -> Tuple[Tuple[str, Path], ...]:
    if not values:
        return DEFAULT_VARIANTS
    result = []
    seen = set()
    for value in values:
        if "=" not in value:
            raise SystemExit(f"invalid --variant {value!r}; expected LABEL=CSV")
        label, raw_path = value.split("=", 1)
        label = label.strip()
        if not label or label == "full" or label in seen:
            raise SystemExit(f"invalid or duplicate variant label: {label!r}")
        seen.add(label)
        result.append((label, Path(raw_path).expanduser()))
    return tuple(result)


def read_csv(path: Path) -> Dict[Tuple[str, str], Dict[str, object]]:
    if not path.is_file():
        raise SystemExit(f"benchmark CSV not found: {path}")
    result: Dict[Tuple[str, str], Dict[str, object]] = {}
    failures = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or ())
        if missing:
            raise SystemExit(f"{path} is missing columns: {sorted(missing)}")
        for line_number, row in enumerate(reader, start=2):
            dataset = row["benchmark"].strip().lower()
            pair_id = row["id"].strip()
            key = (dataset, pair_id)
            if key in result:
                raise SystemExit(f"duplicate pair {key} in {path}:{line_number}")
            if row["status"].strip().lower() != "ok":
                failures.append(key)
                continue
            try:
                result[key] = {
                    metric: float(row[column])
                    for metric, column in METRIC_COLUMNS.items()
                }
            except ValueError as exc:
                raise SystemExit(
                    f"non-numeric metric for {key} in {path}:{line_number}"
                ) from exc
    if failures:
        preview = ", ".join(f"{d}/{i}" for d, i in failures[:5])
        raise SystemExit(
            f"{path} contains {len(failures)} non-ok rows (e.g. {preview}); "
            "paired statistics require a complete result set"
        )
    return result


def dataset_ids(
    rows: Mapping[Tuple[str, str], Mapping[str, object]], dataset: str
) -> Tuple[str, ...]:
    return tuple(sorted(pair_id for ds, pair_id in rows if ds == dataset))


def validate_pairs(
    full: Mapping[Tuple[str, str], Mapping[str, object]],
    variant: Mapping[Tuple[str, str], Mapping[str, object]],
    label: str,
    dataset: str,
) -> Tuple[str, ...]:
    full_ids = dataset_ids(full, dataset)
    variant_ids = dataset_ids(variant, dataset)
    if not full_ids:
        raise SystemExit(f"Full EAlign has no rows for dataset {dataset!r}")
    if full_ids != variant_ids:
        missing = sorted(set(full_ids) - set(variant_ids))
        extra = sorted(set(variant_ids) - set(full_ids))
        raise SystemExit(
            f"pair mismatch for {dataset}/{label}: "
            f"missing={missing[:5]} ({len(missing)} total), "
            f"extra={extra[:5]} ({len(extra)} total)"
        )
    return full_ids


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


def macro_summary(
    rows: Mapping[Tuple[str, str], Mapping[str, object]],
    dataset: str,
    ids: Iterable[str],
    bad_threshold: float,
) -> Dict[str, object]:
    ids = tuple(ids)
    f1 = np.array([rows[(dataset, pair_id)]["f1"] for pair_id in ids])
    return {
        "pairs": len(ids),
        "uppercase_macro_precision": float(
            np.mean([rows[(dataset, i)]["precision"] for i in ids])
        ),
        "uppercase_macro_recall": float(
            np.mean([rows[(dataset, i)]["recall"] for i in ids])
        ),
        "uppercase_macro_f1": float(f1.mean()),
        "uppercase_bad_cases_f1_lt_threshold": int(
            np.count_nonzero(f1 < bad_threshold)
        ),
    }


def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    variants = parse_variants(args.variant)
    datasets = tuple(d.lower() for d in args.dataset) or DEFAULT_DATASETS
    full = read_csv(args.full_csv)
    variant_rows = {label: read_csv(path) for label, path in variants}

    summary_rows = []
    comparison_rows = []
    for dataset in datasets:
        ids = dataset_ids(full, dataset)
        if not ids:
            raise SystemExit(f"no Full EAlign rows for dataset {dataset!r}")
        full_summary = macro_summary(full, dataset, ids, args.bad_threshold)
        summary_rows.append(
            {"dataset": dataset, "variant": "full", **full_summary, "delta_f1_pct_points": 0.0}
        )
        for label, _ in variants:
            current = variant_rows[label]
            paired_ids = validate_pairs(full, current, label, dataset)
            current_summary = macro_summary(
                current, dataset, paired_ids, args.bad_threshold
            )
            delta = (
                current_summary["uppercase_macro_f1"]
                - full_summary["uppercase_macro_f1"]
            ) * 100.0
            summary_rows.append(
                {
                    "dataset": dataset,
                    "variant": label,
                    **current_summary,
                    "delta_f1_pct_points": delta,
                }
            )
            differences = 100.0 * np.array(
                [
                    full[(dataset, pair_id)]["f1"]
                    - current[(dataset, pair_id)]["f1"]
                    for pair_id in paired_ids
                ],
                dtype=float,
            )
            low, high = paired_bootstrap_ci(
                differences,
                args.bootstrap,
                args.confidence,
                args.seed,
            )
            comparison_rows.append(
                {
                    "dataset": dataset,
                    "variant": label,
                    "pairs": len(paired_ids),
                    "variant_minus_full_f1_pct_points": delta,
                    "full_minus_variant_f1_pct_points": float(differences.mean()),
                    "ci_low_pct_points": low,
                    "ci_high_pct_points": high,
                }
            )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.output_dir / "ealign_score_ablation_summary.csv"
    comparisons_path = args.output_dir / "ealign_score_ablation_paired_ci.csv"
    json_path = args.output_dir / "ealign_score_ablation_statistics.json"
    write_csv(summary_path, summary_rows)
    write_csv(comparisons_path, comparison_rows)
    json_path.write_text(
        json.dumps(
            {
                "definitions": {
                    "reference": (
                        "only released alignment columns in which both "
                        "sequences are uppercase"
                    ),
                    "uppercase_macro_metrics": (
                        "arithmetic mean of per-pair uppercase precision, "
                        "recall, or F1; stored as fractions in [0, 1]"
                    ),
                    "delta_f1": (
                        "uppercase macro F1(variant) - uppercase macro "
                        "F1(full), percentage points"
                    ),
                    "paired_difference": (
                        "per-pair uppercase F1(full) - uppercase "
                        "F1(variant), percentage points"
                    ),
                    "bad_case": (
                        f"per-pair uppercase F1 < {args.bad_threshold}"
                    ),
                    "ci": "percentile paired bootstrap CI of the mean paired difference",
                },
                "bootstrap": {
                    "resamples": args.bootstrap,
                    "seed": args.seed,
                    "confidence": args.confidence,
                    "rng": "numpy.random.default_rng (PCG64)",
                },
                "inputs": {
                    "full": str(args.full_csv.resolve()),
                    "variants": {
                        label: str(path.resolve()) for label, path in variants
                    },
                },
                "summary": summary_rows,
                "paired_comparisons": comparison_rows,
            },
            indent=2,
        )
        + "\n"
    )

    print(
        "dataset   variant              P       R      F1     DeltaF1  Bad   "
        "Full-variant  95% CI"
    )
    comparison_index = {
        (row["dataset"], row["variant"]): row for row in comparison_rows
    }
    for row in summary_rows:
        key = (row["dataset"], row["variant"])
        if row["variant"] == "full":
            effect = "--            --"
        else:
            comparison = comparison_index[key]
            effect = (
                f"{comparison['full_minus_variant_f1_pct_points']:7.2f}  "
                f"[{comparison['ci_low_pct_points']:.2f}, "
                f"{comparison['ci_high_pct_points']:.2f}]"
            )
        print(
            f"{row['dataset']:<9} {row['variant']:<17} "
            f"{100.0 * row['uppercase_macro_precision']:7.2f} "
            f"{100.0 * row['uppercase_macro_recall']:7.2f} "
            f"{100.0 * row['uppercase_macro_f1']:7.2f} "
            f"{row['delta_f1_pct_points']:+9.2f} "
            f"{row['uppercase_bad_cases_f1_lt_threshold']:4d}   {effect}"
        )
    print(f"\nWrote {summary_path}")
    print(f"Wrote {comparisons_path}")
    print(f"Wrote {json_path}")


if __name__ == "__main__":
    main()
