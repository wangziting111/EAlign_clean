#!/usr/bin/env python3
"""Build SABmark SupFP benchmark, alignment-job, and query-target indexes."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


EXPECTED = {
    "rows": 79095,
    "TT": 19092,
    "TF": 41116,
    "FF": 18887,
    "benchmark_pairs": 60208,
    "query_targets": 79300,
    "true_queries": 3280,
}


def clean_sid(value: str) -> str:
    sid = value.rsplit(":", 1)[-1]
    if len(sid) != 7 or not sid.startswith("d"):
        raise ValueError(f"unexpected SCOP SID: {value!r}")
    return sid


def parse_meta(value):
    return json.loads(value) if isinstance(value, str) else value


def write_tsv(path: Path, fields, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jsonl", type=Path, default=here / "sup_fp.jsonl")
    parser.add_argument("--structures", type=Path, default=here.parent / "structures")
    parser.add_argument("--output-dir", type=Path, default=here / "manifests")
    args = parser.parse_args()

    jsonl = args.jsonl.resolve(strict=True)
    structures = args.structures.resolve()
    output_dir = args.output_dir.resolve()

    records = []
    type_counts = Counter()
    membership = defaultdict(dict)
    with jsonl.open(encoding="utf-8") as handle:
        for row_number, line in enumerate(handle, 1):
            record = json.loads(line)
            meta = parse_meta(record["meta"])
            sid1 = clean_sid(record["seq1_id"])
            sid2 = clean_sid(record["seq2_id"])
            flag1 = int(meta["true_pos1"])
            flag2 = int(meta["true_pos2"])
            pair_type = "TT" if flag1 and flag2 else "FF" if not flag1 and not flag2 else "TF"
            type_counts[pair_type] += 1
            group = record["group_id"]
            for sid, flag in ((sid1, flag1), (sid2, flag2)):
                previous = membership[group].get(sid)
                if previous is not None and previous != flag:
                    raise RuntimeError(f"inconsistent membership for {group}:{sid}")
                membership[group][sid] = flag
            if pair_type != "FF":
                records.append({
                    "source_row": row_number,
                    "pair_id": record["pair_id"],
                    "group_id": group,
                    "sid1": sid1,
                    "sid2": sid2,
                    "true_pos1": flag1,
                    "true_pos2": flag2,
                    "pair_type": pair_type,
                    "label": 1 if pair_type == "TT" else 0,
                    "percent_identity": record.get("percent_identity", ""),
                    "ref_alignment_length": len(record.get("ref_alignment") or []),
                })

    if sum(type_counts.values()) != EXPECTED["rows"]:
        raise RuntimeError(f"expected {EXPECTED['rows']} rows, got {sum(type_counts.values())}")
    for pair_type in ("TT", "TF", "FF"):
        if type_counts[pair_type] != EXPECTED[pair_type]:
            raise RuntimeError(
                f"expected {EXPECTED[pair_type]} {pair_type}, got {type_counts[pair_type]}"
            )
    if len(records) != EXPECTED["benchmark_pairs"]:
        raise RuntimeError(f"expected {EXPECTED['benchmark_pairs']} TT+TF rows, got {len(records)}")

    canonical_keys = sorted({tuple(sorted((r["sid1"], r["sid2"]))) for r in records})
    job_for_key = {key: f"job{index:06d}" for index, key in enumerate(canonical_keys, 1)}
    occurrences = Counter(tuple(sorted((r["sid1"], r["sid2"]))) for r in records)

    benchmark_rows = []
    query_rows = []
    for pair_index, record in enumerate(records, 1):
        key = tuple(sorted((record["sid1"], record["sid2"])))
        job_id = job_for_key[key]
        row = dict(record)
        row.update({
            "pair_index": pair_index,
            "job_id": job_id,
            "canonical_sid_a": key[0],
            "canonical_sid_b": key[1],
            "cross_group_occurrences": occurrences[key],
        })
        benchmark_rows.append(row)

        if record["pair_type"] == "TT":
            directions = (
                (record["sid1"], record["sid2"]),
                (record["sid2"], record["sid1"]),
            )
        elif record["true_pos1"]:
            directions = ((record["sid1"], record["sid2"]),)
        else:
            directions = ((record["sid2"], record["sid1"]),)
        for query, target in directions:
            query_rows.append({
                "query_index": len(query_rows) + 1,
                "group_id": record["group_id"],
                "query_sid": query,
                "target_sid": target,
                "label": record["label"],
                "pair_type": record["pair_type"],
                "pair_index": pair_index,
                "job_id": job_id,
            })

    job_rows = []
    for sid_a, sid_b in canonical_keys:
        key = (sid_a, sid_b)
        job_rows.append({
            "job_id": job_for_key[key],
            "sid_a": sid_a,
            "sid_b": sid_b,
            "pdb_a": str(structures / sid_a[2:4] / f"{sid_a}.ent"),
            "pdb_b": str(structures / sid_b[2:4] / f"{sid_b}.ent"),
            "benchmark_occurrences": occurrences[key],
        })

    true_queries = {row["query_sid"] for row in query_rows}
    if len(query_rows) != EXPECTED["query_targets"]:
        raise RuntimeError(f"expected {EXPECTED['query_targets']} query targets, got {len(query_rows)}")
    if len(true_queries) != EXPECTED["true_queries"]:
        raise RuntimeError(f"expected {EXPECTED['true_queries']} queries, got {len(true_queries)}")

    pair_fields = (
        "pair_index", "job_id", "pair_id", "group_id", "sid1", "sid2",
        "true_pos1", "true_pos2", "pair_type", "label", "canonical_sid_a",
        "canonical_sid_b", "cross_group_occurrences", "percent_identity",
        "ref_alignment_length", "source_row",
    )
    job_fields = ("job_id", "sid_a", "sid_b", "pdb_a", "pdb_b", "benchmark_occurrences")
    query_fields = (
        "query_index", "group_id", "query_sid", "target_sid", "label",
        "pair_type", "pair_index", "job_id",
    )
    write_tsv(output_dir / "supfp_tt_tf_pairs.tsv", pair_fields, benchmark_rows)
    write_tsv(output_dir / "supfp_ealign_jobs.tsv", job_fields, job_rows)
    write_tsv(output_dir / "supfp_query_targets.tsv", query_fields, query_rows)

    missing = sorted({
        sid
        for row in job_rows
        for sid, field in ((row["sid_a"], "pdb_a"), (row["sid_b"], "pdb_b"))
        if not Path(row[field]).is_file()
    })
    summary = {
        "source_rows": sum(type_counts.values()),
        "pair_type_counts": dict(sorted(type_counts.items())),
        "benchmark_tt_tf_rows": len(benchmark_rows),
        "globally_unique_ealign_jobs": len(job_rows),
        "cross_group_duplicate_rows": len(benchmark_rows) - len(job_rows),
        "directed_query_target_rows": len(query_rows),
        "true_queries": len(true_queries),
        "missing_structure_files": len(missing),
    }
    (output_dir / "supfp_pair_index_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    (output_dir / "missing_structure_sids.txt").write_text(
        "".join(f"{sid}\n" for sid in missing), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
