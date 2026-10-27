#!/usr/bin/env python3
"""Run EAlign once for every TT/TF record in a SABmark JSONL file."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


OUTPUT_SUFFIXES = (".fasta", ".score", ".local", ".pdb")
JSONL_NAMES = {"supfp": "sup_fp.jsonl", "twifp": "twi_fp.jsonl"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_records(jsonl_path: Path, subset: str, structures: Path) -> list[dict]:
    """Keep TT/TF rows in source order and preserve seq1 -> seq2 orientation."""
    records = []
    counts = Counter()

    with jsonl_path.open(encoding="utf-8") as handle:
        for source_line, line in enumerate(handle, 1):
            source = json.loads(line)
            meta = json.loads(source["meta"]) if isinstance(source["meta"], str) else source["meta"]
            flags = (int(meta["true_pos1"]), int(meta["true_pos2"]))

            if flags == (0, 0):
                counts["FF"] += 1
                continue
            if flags == (1, 1):
                pair_type, label = "TT", 1
            elif flags == (1, 0):
                pair_type, label = "TF", 0
            else:
                raise ValueError(
                    f"{jsonl_path}:{source_line}: unsupported truth flags {flags}"
                )
            counts[pair_type] += 1

            sid1 = source["seq1_id"].split(":", 1)[-1]
            sid2 = source["seq2_id"].split(":", 1)[-1]
            records.append(
                {
                    "record_id": f"{subset}_line{source_line:06d}",
                    "source_line": source_line,
                    "pair_id": source["pair_id"],
                    "group_id": source["group_id"],
                    "pair_type": pair_type,
                    "label": label,
                    "sid1": sid1,
                    "sid2": sid2,
                    "pdb1": str((structures / sid1[2:4] / f"{sid1}.ent").resolve()),
                    "pdb2": str((structures / sid2[2:4] / f"{sid2}.ent").resolve()),
                }
            )

    print(
        f"input={jsonl_path} TT={counts['TT']} TF={counts['TF']} "
        f"FF_skipped={counts['FF']} selected={len(records)}",
        flush=True,
    )
    return records


def parse_score(path: Path) -> dict:
    for line in path.read_text(errors="replace").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 21 or fields[4] != "->" or fields[8] != "->":
            raise ValueError(f"unexpected EAlign score record: {line}")
        return {
            "name1": fields[0],
            "name2": fields[1],
            "len1": int(fields[2]),
            "len2": int(fields[3]),
            "blosum": float(fields[5]),
            "score_3di": float(fields[6]),
            "combined_score": float(fields[7]),
            "lali": int(fields[9]),
            "rmsd": float(fields[10]),
            "tmscore": float(fields[11]),
            "maxsub": float(fields[13]),
            "gdt_ts": float(fields[14]),
            "gdt_ha": float(fields[15]),
            "seqid": float(fields[17]),
            "normalization_length": int(fields[18]),
            "dcutoff": float(fields[19]),
            "ugdt": float(fields[20]),
        }
    raise ValueError(f"no EAlign score record in {path}")


def record_directory(output_root: Path, record: dict) -> Path:
    bucket = f"{record['source_line'] // 1000:03d}"
    return output_root / "records" / bucket / record["record_id"]


def write_json_atomic(path: Path, value: dict) -> None:
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def reusable(result: dict, destination: Path, record: dict, binary_hash: str) -> bool:
    return (
        result.get("status") == "ok"
        and result.get("binary_sha256") == binary_hash
        and result.get("pdb1") == record["pdb1"]
        and result.get("pdb2") == record["pdb2"]
        and all((destination / f"{record['record_id']}{suffix}").is_file()
                for suffix in OUTPUT_SUFFIXES)
    )


def run_one(
    record: dict,
    binary: Path,
    binary_hash: str,
    output_root: Path,
    overwrite: bool,
) -> dict:
    destination = record_directory(output_root, record)
    destination.mkdir(parents=True, exist_ok=True)
    result_path = destination / "result.json"

    if result_path.is_file() and not overwrite:
        previous = json.loads(result_path.read_text())
        if reusable(previous, destination, record, binary_hash):
            previous["cached"] = True
            return previous

    result = {
        **record,
        "binary": str(binary),
        "binary_sha256": binary_hash,
        "status": "error",
        "cached": False,
        "error": "",
    }

    for field in ("pdb1", "pdb2"):
        if not Path(record[field]).is_file():
            result["error"] = f"missing input: {record[field]}"
            write_json_atomic(result_path, result)
            return result

    record_id = record["record_id"]
    command = [str(binary), record["pdb1"], record["pdb2"], "-p", "1", "-o", record_id]
    result["command"] = command

    with tempfile.TemporaryDirectory(prefix=f"ealign_{record_id}_") as scratch_name:
        scratch = Path(scratch_name)
        started = time.perf_counter()
        try:
            completed = subprocess.run(
                command,
                cwd=scratch,
                capture_output=True,
                text=True,
            )
            alignment_seconds = time.perf_counter() - started
            result["alignment_seconds"] = round(alignment_seconds, 6)
            result["returncode"] = completed.returncode
            (destination / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
            (destination / "stderr.txt").write_text(completed.stderr, encoding="utf-8")

            produced = [scratch / f"{record_id}{suffix}" for suffix in OUTPUT_SUFFIXES]
            if completed.returncode != 0:
                result["error"] = f"returncode={completed.returncode}"
            elif not all(path.is_file() for path in produced):
                result["error"] = "missing one or more expected EAlign output files"
            else:
                for source in produced:
                    shutil.move(str(source), destination / source.name)
                result["metrics"] = parse_score(destination / f"{record_id}.score")
                result["status"] = "ok"
        except Exception as exc:
            result["error"] = f"{type(exc).__name__}: {exc}"

    write_json_atomic(result_path, result)
    return result


def write_summary(output_root: Path, results: list[dict]) -> Path:
    fields = (
        "record_id", "source_line", "pair_id", "group_id", "pair_type", "label",
        "sid1", "sid2", "pdb1", "pdb2", "status", "cached",
        "alignment_seconds",
        "len1", "len2", "combined_score", "lali", "rmsd", "tmscore",
        "maxsub", "gdt_ts", "gdt_ha", "seqid", "normalization_length",
        "dcutoff", "ugdt", "error", "result_json",
    )
    path = output_root / "ealign_results.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for result in sorted(results, key=lambda row: row["source_line"]):
            metrics = result.get("metrics", {})
            row = {field: result.get(field, metrics.get(field, "")) for field in fields}
            row["result_json"] = str(
                record_directory(output_root, result) / "result.json"
            )
            writer.writerow(row)
    return path


def main() -> int:
    here = Path(__file__).resolve().parent
    repo = here.parents[1]

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("subset", choices=("supfp", "twifp"))
    parser.add_argument("--input-jsonl", type=Path)
    parser.add_argument("--structures", type=Path, default=here / "structures")
    parser.add_argument("--binary", type=Path, default=repo / "binary/ealign")
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    if args.workers < 1:
        parser.error("--workers must be positive")

    jsonl_path = (
        args.input_jsonl or here / args.subset / JSONL_NAMES[args.subset]
    ).resolve(strict=True)
    structures = args.structures.resolve(strict=True)
    binary = args.binary.resolve(strict=True)
    output_root = (
        args.output_root or here / "confirm" / "results" / args.subset / "ealign"
    ).resolve()

    records = load_records(jsonl_path, args.subset, structures)
    binary_hash = sha256(binary)
    print(
        f"method=ealign subset={args.subset} records={len(records)} "
        f"workers={args.workers} output={output_root} binary_sha256={binary_hash}",
        flush=True,
    )

    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                run_one, record, binary, binary_hash, output_root, args.overwrite
            ): record
            for record in records
        }
        for completed_count, future in enumerate(as_completed(futures), 1):
            result = future.result()
            results.append(result)
            if completed_count % 100 == 0 or completed_count == len(records):
                ok = sum(row["status"] == "ok" for row in results)
                cached = sum(bool(row.get("cached")) for row in results)
                print(
                    f"progress={completed_count}/{len(records)} "
                    f"ok={ok} cached={cached}",
                    flush=True,
                )

    summary = write_summary(output_root, results)
    failed = [row for row in results if row["status"] != "ok"]
    print(f"summary={summary}")
    print(f"ok={len(results) - len(failed)} failed={len(failed)}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
