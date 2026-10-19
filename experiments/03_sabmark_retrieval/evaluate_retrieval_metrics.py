#!/usr/bin/env python3
"""Calculate query-macro AUROC, mAP, and ROC1 for TwiFP or SupFP.

No third-party Python packages are required.

Definitions:
  AUROC: positive-vs-negative pair comparisons; a score tie counts as 0.5.
  AP: area under the precision-recall step curve, grouping equal scores.
  ROC1: fraction of positives scoring strictly above the highest negative.
  minlen calibration: raw score / min(length_a, length_b).
  meanlen calibration: raw score / ((length_a + length_b) / 2).
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


METHODS = (
    ("TMalign_minlen", "results/tmalign/tmalign_results.tsv", "tm_score_minlen", "none"),
    # ("TMalign_mean", "results/tmalign/tmalign_results.tsv", "tm_score_mean", "none"),
    # ("DeepAlign_DeepScore_raw", "results/deepalign/deepalign_results.tsv", "deep_score", "none"),
    ("DeepAlign_DeepScore_per_minlen", "results/deepalign/deepalign_results.tsv", "deep_score", "minlen"),
    # ("DeepAlign_DeepScore_per_meanlen", "results/deepalign/deepalign_results.tsv", "deep_score", "meanlen"),
    # ("EAlign_CombinedScore_raw", "results/ealign/ealign_results.tsv", "combined_score", "none"),
    ("EAlign_CombinedScore_per_minlen", "results/ealign/ealign_results.tsv", "combined_score", "minlen"),
    # ("EAlign_CombinedScore_per_meanlen", "results/ealign/ealign_results.tsv", "combined_score", "meanlen"),
)


def read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_lengths(path: Path):
    lengths = {}
    for row in read_tsv(path):
        if row["status"] == "ok":
            lengths[row["job_id"]] = (int(row["length_a"]), int(row["length_b"]))
    return lengths


def load_scores(path: Path, score_column: str, normalization: str, lengths):
    scores = {}
    for row in read_tsv(path):
        if row["status"] != "ok":
            continue
        job_id = row["job_id"]
        score = float(row[score_column])
        if normalization != "none":
            if job_id not in lengths:
                raise RuntimeError(f"missing lengths for {job_id}")
            length_a, length_b = lengths[job_id]
            if normalization == "minlen":
                score /= min(length_a, length_b)
            elif normalization == "meanlen":
                score /= (length_a + length_b) / 2.0
            else:
                raise ValueError(f"unknown normalization: {normalization}")
        scores[job_id] = score
    return scores


def auc(labels, scores):
    positives = [score for label, score in zip(labels, scores) if label == 1]
    negatives = [score for label, score in zip(labels, scores) if label == 0]
    correct = sum(
        (positive > negative) + 0.5 * (positive == negative)
        for positive in positives
        for negative in negatives
    )
    return correct / (len(positives) * len(negatives))


def average_precision(labels, scores):
    # Equal scores enter the ranking together, so the result does not depend
    # on an arbitrary row order inside a tie.
    score_groups = defaultdict(lambda: [0, 0])
    for label, score in zip(labels, scores):
        score_groups[score][0 if label == 1 else 1] += 1

    total_positives = sum(labels)
    true_positives = 0
    false_positives = 0
    result = 0.0
    for score in sorted(score_groups, reverse=True):
        added_positives, added_negatives = score_groups[score]
        true_positives += added_positives
        false_positives += added_negatives
        precision = true_positives / (true_positives + false_positives)
        recall_increase = added_positives / total_positives
        result += recall_increase * precision
    return result


def roc1(labels, scores):
    positives = [score for label, score in zip(labels, scores) if label == 1]
    negatives = [score for label, score in zip(labels, scores) if label == 0]
    highest_negative = max(negatives)
    return sum(score > highest_negative for score in positives) / len(positives)


def evaluate(query_rows, scores):
    queries = defaultdict(lambda: {"labels": [], "scores": []})
    for row in query_rows:
        job_id = row["job_id"]
        if job_id not in scores:
            raise RuntimeError(f"missing score for {job_id}")
        query_key = (row["group_id"], row["query_sid"])
        query = queries[query_key]
        query["labels"].append(int(row["label"]))
        query["scores"].append(scores[job_id])

    values = {"auroc": [], "map": [], "roc1": []}
    excluded = []
    for query_key, data in queries.items():
        labels = data["labels"]
        query_scores = data["scores"]
        if 0 not in labels or 1 not in labels:
            excluded.append(query_key)
            continue
        values["auroc"].append(auc(labels, query_scores))
        values["map"].append(average_precision(labels, query_scores))
        values["roc1"].append(roc1(labels, query_scores))

    count = len(values["auroc"])
    if count == 0:
        raise RuntimeError("no query has both positive and negative targets")
    return {
        "queries": count,
        "total_queries": len(queries),
        "excluded_queries": len(excluded),
        "auroc": sum(values["auroc"]) / count,
        "map": sum(values["map"]) / count,
        "roc1": sum(values["roc1"]) / count,
    }


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    sabmark_root = script_dir.parents[1] / "dataset/SABmark"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-subset", "--subset", choices=("twifp", "supfp"), default="twifp",
        help="SABmark FP subset to evaluate (default: twifp)",
    )
    parser.add_argument(
        "--root", type=Path,
        help="override the subset directory selected by -subset",
    )
    args = parser.parse_args()

    root = (args.root or sabmark_root / args.subset).resolve()
    query_index = root / f"manifests/{args.subset}_query_targets.tsv"
    query_rows = read_tsv(query_index)
    lengths = load_lengths(root / "results/deepalign/deepalign_results.tsv")

    protocol = None
    print("method\tscore\tqueries\tAUROC\tmAP\tROC1")
    for method, relative_path, score_column, normalization in METHODS:
        result_path = root / relative_path
        scores = load_scores(result_path, score_column, normalization, lengths)
        metrics = evaluate(query_rows, scores)
        current_protocol = (
            metrics["total_queries"], metrics["queries"], metrics["excluded_queries"]
        )
        if protocol is None:
            protocol = current_protocol
        elif current_protocol != protocol:
            raise RuntimeError("methods were evaluated on different query sets")
        score_name = score_column if normalization == "none" else f"{score_column}/{normalization}"
        print(
            f"{method}\t{score_name}\t{metrics['queries']}\t"
            f"{metrics['auroc']:.6f}\t{metrics['map']:.6f}\t"
            f"{metrics['roc1']:.6f}"
        )
    total, evaluated, excluded = protocol
    print(
        f"# subset={args.subset} total_queries={total} "
        f"evaluated_queries={evaluated} excluded_single_class={excluded}"
    )


if __name__ == "__main__":
    main()
