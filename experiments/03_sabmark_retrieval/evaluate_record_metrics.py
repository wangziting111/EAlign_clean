#!/usr/bin/env python3
"""Evaluate record-level SABmark TT/TF results with query-macro metrics.

The default paths follow the clean release layout.  Command-line overrides are
kept so the evaluator can also be run against an independently prepared copy.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


JSONL_NAMES = {"supfp": "sup_fp.jsonl", "twifp": "twi_fp.jsonl"}
EXPECTED_RECORDS = {"supfp": 60208, "twifp": 33518}
EXPECTED_RELATIONS = {"supfp": 79300, "twifp": 44185}

METHODS = (
    {
        "name": "TMalign_minlen",
        "directory": "tmalign",
        "filename": "tmalign_results.tsv",
        "score_column": "tm_score_minlen",
        "length_columns": None,
        "score_definition": "tm_score_minlen",
    },
    {
        "name": "DeepAlign_DeepScore_per_minlen",
        "directory": "deepalign",
        "filename": "deepalign_results.tsv",
        "score_column": "deep_score",
        "length_columns": ("length1", "length2"),
        "score_definition": "deep_score/min(length1,length2)",
    },
    {
        "name": "EAlign_CombinedScore_per_minlen",
        "directory": "ealign",
        "filename": "ealign_results.tsv",
        "score_column": "combined_score",
        "length_columns": ("len1", "len2"),
        "score_definition": "combined_score/min(len1,len2)",
    },
    {
        "name": "GTalign_CPU_TMscore_minlen",
        "directory": "gtalign_cpu_patched_t2_w16",
        "filename": "gtalign_results.tsv",
        "score_column": "tm_score_minlen",
        "length_columns": None,
        "score_definition": "tm_score_minlen",
    },
    {
        "name": "Foldseek_bit_score",
        "directory": "foldseek",
        "filename": "foldseek_results.tsv",
        "score_column": "bits",
        "length_columns": None,
        "score_definition": "bits",
    },
    {
        "name": "Foldseek_probability",
        "directory": "foldseek",
        "filename": "foldseek_results.tsv",
        "score_column": "prob",
        "length_columns": None,
        "score_definition": "prob",
    },
)


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def load_records(jsonl_path: Path, subset: str) -> dict[str, dict]:
    """Reconstruct retained records directly from JSONL source lines."""
    records = {}

    with jsonl_path.open(encoding="utf-8") as handle:
        for source_line, line in enumerate(handle, 1):
            source = json.loads(line)
            meta = json.loads(source["meta"]) if isinstance(source["meta"], str) else source["meta"]
            flags = (int(meta["true_pos1"]), int(meta["true_pos2"]))

            if flags == (0, 0):
                continue
            if flags == (1, 1):
                pair_type, label = "TT", 1
            elif flags == (1, 0):
                pair_type, label = "TF", 0
            else:
                raise ValueError(
                    f"{jsonl_path}:{source_line}: unsupported truth flags {flags}"
                )

            record_id = f"{subset}_line{source_line:06d}"
            sid1 = source["seq1_id"].split(":", 1)[-1]
            sid2 = source["seq2_id"].split(":", 1)[-1]
            records[record_id] = {
                "record_id": record_id,
                "source_line": source_line,
                "pair_id": source["pair_id"],
                "group_id": source["group_id"],
                "pair_type": pair_type,
                "label": label,
                "sid1": sid1,
                "sid2": sid2,
            }

    expected = EXPECTED_RECORDS[subset]
    if len(records) != expected:
        raise RuntimeError(
            f"{subset}: retained {len(records)} records, expected {expected}"
        )
    return records


def build_relations(records: dict[str, dict]) -> list[dict]:
    relations = []

    for record in records.values():
        common = {
            "record_id": record["record_id"],
            "group_id": record["group_id"],
        }
        if record["pair_type"] == "TT":
            relations.append(
                {
                    **common,
                    "query_sid": record["sid1"],
                    "target_sid": record["sid2"],
                    "label": 1,
                }
            )
            relations.append(
                {
                    **common,
                    "query_sid": record["sid2"],
                    "target_sid": record["sid1"],
                    "label": 1,
                }
            )
        else:
            relations.append(
                {
                    **common,
                    "query_sid": record["sid1"],
                    "target_sid": record["sid2"],
                    "label": 0,
                }
            )

    return relations


def load_method_scores(
    path: Path,
    method: dict,
    records: dict[str, dict],
) -> dict[str, float]:
    rows = read_tsv(path)
    by_id = {}

    for row_number, row in enumerate(rows, 2):
        record_id = row["record_id"]
        if record_id in by_id:
            raise RuntimeError(f"{path}:{row_number}: duplicate record_id {record_id}")
        by_id[record_id] = row

    expected_ids = set(records)
    observed_ids = set(by_id)
    missing = expected_ids - observed_ids
    extra = observed_ids - expected_ids
    if missing or extra:
        raise RuntimeError(
            f"{path}: record_id mismatch: missing={len(missing)} extra={len(extra)}"
        )

    scores = {}
    for record_id, record in records.items():
        row = by_id[record_id]
        if row["status"] != "ok":
            raise RuntimeError(f"{path}: {record_id} status={row['status']}")

        for column, expected in (
            ("group_id", record["group_id"]),
            ("pair_type", record["pair_type"]),
            ("sid1", record["sid1"]),
            ("sid2", record["sid2"]),
        ):
            if row[column] != str(expected):
                raise RuntimeError(
                    f"{path}: {record_id} {column}={row[column]!r}, "
                    f"expected {expected!r}"
                )

        score = float(row[method["score_column"]])
        if method["length_columns"] is not None:
            left, right = method["length_columns"]
            length1, length2 = int(row[left]), int(row[right])
            if length1 <= 0 or length2 <= 0:
                raise RuntimeError(f"{path}: invalid lengths for {record_id}")
            score /= min(length1, length2)

        if not math.isfinite(score):
            raise RuntimeError(f"{path}: non-finite score for {record_id}")
        scores[record_id] = score

    return scores


def auc(labels: list[int], scores: list[float]) -> float:
    positives = [score for label, score in zip(labels, scores) if label == 1]
    negatives = [score for label, score in zip(labels, scores) if label == 0]
    correct = sum(
        (positive > negative) + 0.5 * (positive == negative)
        for positive in positives
        for negative in negatives
    )
    return correct / (len(positives) * len(negatives))


def average_precision(labels: list[int], scores: list[float]) -> float:
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
        result += (added_positives / total_positives) * precision

    return result


def roc1(labels: list[int], scores: list[float]) -> float:
    positives = [score for label, score in zip(labels, scores) if label == 1]
    negatives = [score for label, score in zip(labels, scores) if label == 0]
    highest_negative = max(negatives)
    return sum(score > highest_negative for score in positives) / len(positives)


def evaluate_method(
    subset: str,
    method: dict,
    relations: list[dict],
    scores: dict[str, float],
) -> tuple[dict, list[dict]]:
    queries = defaultdict(lambda: {"labels": [], "scores": [], "targets": []})

    for relation in relations:
        query_key = (relation["group_id"], relation["query_sid"])
        query = queries[query_key]
        query["labels"].append(relation["label"])
        query["scores"].append(scores[relation["record_id"]])
        query["targets"].append(relation["target_sid"])

    metric_values = {"auroc": [], "average_precision": [], "roc1": []}
    query_rows = []

    for (group_id, query_sid), data in sorted(queries.items()):
        labels = data["labels"]
        query_scores = data["scores"]
        positives = sum(labels)
        negatives = len(labels) - positives
        excluded = positives == 0 or negatives == 0

        row = {
            "subset": subset,
            "method": method["name"],
            "group_id": group_id,
            "query_sid": query_sid,
            "positives": positives,
            "negatives": negatives,
            "excluded_single_class": int(excluded),
            "auroc": "",
            "average_precision": "",
            "roc1": "",
        }

        if not excluded:
            values = {
                "auroc": auc(labels, query_scores),
                "average_precision": average_precision(labels, query_scores),
                "roc1": roc1(labels, query_scores),
            }
            for key, value in values.items():
                metric_values[key].append(value)
                row[key] = f"{value:.12f}"

        query_rows.append(row)

    evaluated = len(metric_values["auroc"])
    if evaluated == 0:
        raise RuntimeError(f"{subset}/{method['name']}: no evaluable query")

    summary = {
        "subset": subset,
        "method": method["name"],
        "score_definition": method["score_definition"],
        "alignment_records": len(scores),
        "query_target_relations": len(relations),
        "total_queries": len(queries),
        "evaluated_queries": evaluated,
        "excluded_queries": len(queries) - evaluated,
        "auroc": f"{sum(metric_values['auroc']) / evaluated:.12f}",
        "map": f"{sum(metric_values['average_precision']) / evaluated:.12f}",
        "roc1": f"{sum(metric_values['roc1']) / evaluated:.12f}",
    }
    return summary, query_rows


def main() -> None:
    script = Path(__file__).resolve()
    repo_root = script.parents[2]
    sabmark_root = repo_root / "dataset" / "SABmark"
    experiment_result_root = (
        repo_root / "test_results" / "03_sabmark_retrieval"
    )

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("subset", choices=("supfp", "twifp"))
    parser.add_argument("--input-jsonl", type=Path)
    parser.add_argument("--results-root", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    jsonl_path = (
        args.input_jsonl
        or sabmark_root / args.subset / JSONL_NAMES[args.subset]
    ).resolve(strict=True)
    results_root = (
        args.results_root
        or experiment_result_root / "results" / args.subset
    ).resolve(strict=True)
    output_dir = (
        args.output_dir
        or experiment_result_root / "evaluation"
    ).resolve()

    records = load_records(jsonl_path, args.subset)
    relations = build_relations(records)
    if len(relations) != EXPECTED_RELATIONS[args.subset]:
        raise RuntimeError(
            f"{args.subset}: built {len(relations)} relations, "
            f"expected {EXPECTED_RELATIONS[args.subset]}"
        )

    summary_rows = []
    all_query_rows = []

    for method in METHODS:
        result_path = (
            results_root / method["directory"] / method["filename"]
        ).resolve(strict=True)
        scores = load_method_scores(result_path, method, records)
        summary, query_rows = evaluate_method(
            args.subset, method, relations, scores
        )
        summary_rows.append(summary)
        all_query_rows.extend(query_rows)

    metrics_path = output_dir / f"{args.subset}_metrics.tsv"
    query_path = output_dir / f"{args.subset}_query_metrics.tsv"

    write_tsv(
        metrics_path,
        [
            "subset", "method", "score_definition", "alignment_records",
            "query_target_relations", "total_queries", "evaluated_queries",
            "excluded_queries", "auroc", "map", "roc1",
        ],
        summary_rows,
    )
    write_tsv(
        query_path,
        [
            "subset", "method", "group_id", "query_sid", "positives",
            "negatives", "excluded_single_class", "auroc",
            "average_precision", "roc1",
        ],
        all_query_rows,
    )

    print("method\tscore\trecords\trelations\tqueries\tAUROC\tmAP\tROC1")
    for row in summary_rows:
        print(
            f"{row['method']}\t{row['score_definition']}\t"
            f"{row['alignment_records']}\t{row['query_target_relations']}\t"
            f"{row['evaluated_queries']}\t{float(row['auroc']):.6f}\t"
            f"{float(row['map']):.6f}\t{float(row['roc1']):.6f}"
        )
    print(f"metrics={metrics_path}")
    print(f"query_metrics={query_path}")


if __name__ == "__main__":
    main()
