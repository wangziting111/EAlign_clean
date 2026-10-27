#!/usr/bin/env python3
"""Run Foldseek once for every TT/TF record in a SABmark JSONL file."""

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


JSONL_NAMES = {"supfp": "sup_fp.jsonl", "twifp": "twi_fp.jsonl"}
FORMAT_FIELDS = (
    "query",
    "target",
    "bits",
    "evalue",
    "prob",
    "alnlen",
    "qstart",
    "qend",
    "tstart",
    "tend",
    "qlen",
    "tlen",
)
FORMAT_OUTPUT = ",".join(FORMAT_FIELDS)
NORMALIZATION = "first_model_first_protein_chain_column_66_v2"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def foldseek_version(binary: Path) -> str:
    completed = subprocess.run(
        [str(binary), "version"],
        capture_output=True,
        text=True,
        check=True,
    )
    return completed.stdout.strip()


def load_records(jsonl_path: Path, subset: str, structures: Path) -> list[dict]:
    """Keep TT/TF rows in source order and preserve seq1 -> seq2 orientation."""
    records = []
    counts = Counter()

    with jsonl_path.open(encoding="utf-8") as handle:
        for source_line, line in enumerate(handle, 1):
            source = json.loads(line)
            meta = (
                json.loads(source["meta"])
                if isinstance(source["meta"], str)
                else source["meta"]
            )
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
                    "pdb1": str(
                        (structures / sid1[2:4] / f"{sid1}.ent").resolve()
                    ),
                    "pdb2": str(
                        (structures / sid2[2:4] / f"{sid2}.ent").resolve()
                    ),
                }
            )

    print(
        f"input={jsonl_path} TT={counts['TT']} TF={counts['TF']} "
        f"FF_skipped={counts['FF']} selected={len(records)}",
        flush=True,
    )
    return records


def write_foldseek_pdb(source: Path, destination: Path) -> None:
    """Write the first protein chain from the first model as a clean PDB."""
    lines = source.read_bytes().splitlines()
    has_models = any(line[:6].strip() == b"MODEL" for line in lines)

    model_lines = []
    inside_first_model = not has_models
    first_model_seen = False
    for line in lines:
        record_name = line[:6].strip()
        if has_models and record_name == b"MODEL":
            if first_model_seen:
                break
            first_model_seen = True
            inside_first_model = True
            continue
        if has_models and inside_first_model and record_name == b"ENDMDL":
            break
        if inside_first_model:
            model_lines.append(line)

    first_chain = None
    for line in model_lines:
        if (
            line.startswith(b"ATOM  ")
            and len(line) >= 22
            and line[12:16].strip() == b"CA"
        ):
            first_chain = line[21:22]
            break
    if first_chain is None:
        raise ValueError(f"no protein CA atom in first model: {source}")

    coordinate_lines = []
    chain_started = False
    for line in model_lines:
        record_name = line[:6].strip()
        if record_name == b"TER" and chain_started:
            break
        if not line.startswith((b"ATOM  ", b"HETATM")):
            continue
        if len(line) < 22 or line[21:22] != first_chain:
            if chain_started and line.startswith(b"ATOM  "):
                break
            continue
        chain_started = True
        coordinate_lines.append(line[:66])

    if not coordinate_lines:
        raise ValueError(f"no coordinates for first protein chain: {source}")

    with destination.open("wb") as output_handle:
        for line in coordinate_lines:
            output_handle.write(line + b"\n")
        output_handle.write(b"TER\nEND\n")


def parse_result(path: Path) -> dict:
    lines = [line for line in path.read_text().splitlines() if line.strip()]
    if len(lines) != 1:
        raise ValueError(
            f"expected exactly one Foldseek result row, found {len(lines)}"
        )

    values = lines[0].split("\t")
    if len(values) != len(FORMAT_FIELDS):
        raise ValueError(
            f"expected {len(FORMAT_FIELDS)} Foldseek fields, found {len(values)}"
        )

    row = dict(zip(FORMAT_FIELDS, values))
    return {
        "query": row["query"],
        "target": row["target"],
        "bits": float(row["bits"]),
        "evalue": float(row["evalue"]),
        "prob": float(row["prob"]),
        "alnlen": int(row["alnlen"]),
        "qstart": int(row["qstart"]),
        "qend": int(row["qend"]),
        "tstart": int(row["tstart"]),
        "tend": int(row["tend"]),
        "qlen": int(row["qlen"]),
        "tlen": int(row["tlen"]),
    }


def record_directory(output_root: Path, record: dict) -> Path:
    bucket = f"{record['source_line'] // 1000:03d}"
    return output_root / "records" / bucket / record["record_id"]


def write_json_atomic(path: Path, value: dict) -> None:
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def reusable(
    result: dict,
    destination: Path,
    record: dict,
    binary_hash: str,
) -> bool:
    return (
        result.get("status") == "ok"
        and result.get("binary_sha256") == binary_hash
        and result.get("normalization") == NORMALIZATION
        and result.get("pdb1") == record["pdb1"]
        and result.get("pdb2") == record["pdb2"]
        and (destination / "foldseek.tsv").is_file()
        and bool(result.get("metrics"))
    )


def run_one(
    record: dict,
    binary: Path,
    binary_hash: str,
    version: str,
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
        "foldseek_version": version,
        "normalization": NORMALIZATION,
        "status": "error",
        "cached": False,
        "error": "",
    }

    for field in ("pdb1", "pdb2"):
        if not Path(record[field]).is_file():
            result["error"] = f"missing input: {record[field]}"
            write_json_atomic(result_path, result)
            return result

    try:
        with tempfile.TemporaryDirectory(
            prefix=f"foldseek_{record['record_id']}_"
        ) as scratch_name:
            scratch = Path(scratch_name)
            query_pdb = scratch / f"{record['sid1']}.pdb"
            target_pdb = scratch / f"{record['sid2']}.pdb"
            output_tsv = scratch / "foldseek.tsv"
            temporary_directory = scratch / "tmp"

            write_foldseek_pdb(Path(record["pdb1"]), query_pdb)
            write_foldseek_pdb(Path(record["pdb2"]), target_pdb)

            command = [
                str(binary),
                "easy-search",
                str(query_pdb),
                str(target_pdb),
                str(output_tsv),
                str(temporary_directory),
                "--alignment-type",
                "2",
                "--exhaustive-search",
                "1",
                "-e",
                "inf",
                "--threads",
                "1",
                "--format-output",
                FORMAT_OUTPUT,
            ]
            result["command"] = command

            started = time.perf_counter()
            completed = subprocess.run(
                command,
                cwd=scratch,
                capture_output=True,
                text=True,
            )
            alignment_seconds = time.perf_counter() - started

            result["alignment_seconds"] = round(alignment_seconds, 6)
            result["returncode"] = completed.returncode
            (destination / "stdout.txt").write_text(
                completed.stdout, encoding="utf-8"
            )
            (destination / "stderr.txt").write_text(
                completed.stderr, encoding="utf-8"
            )

            if completed.returncode != 0:
                result["error"] = f"returncode={completed.returncode}"
            elif not output_tsv.is_file():
                result["error"] = "missing Foldseek output TSV"
            else:
                result["metrics"] = parse_result(output_tsv)
                shutil.copy2(output_tsv, destination / "foldseek.tsv")
                result["status"] = "ok"
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"

    write_json_atomic(result_path, result)
    return result


def write_summary(output_root: Path, results: list[dict]) -> Path:
    fields = (
        "record_id",
        "source_line",
        "pair_id",
        "group_id",
        "pair_type",
        "label",
        "sid1",
        "sid2",
        "pdb1",
        "pdb2",
        "status",
        "cached",
        "alignment_seconds",
        "query",
        "target",
        "bits",
        "evalue",
        "prob",
        "alnlen",
        "qstart",
        "qend",
        "tstart",
        "tend",
        "qlen",
        "tlen",
        "normalization",
        "error",
        "result_json",
    )
    path = output_root / "foldseek_results.tsv"
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
    structure_alignment_root = here.parents[2]

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("subset", choices=("supfp", "twifp"))
    parser.add_argument("--input-jsonl", type=Path)
    parser.add_argument("--structures", type=Path, default=here / "structures")
    parser.add_argument(
        "--binary",
        type=Path,
        default=(
            structure_alignment_root
            / "foldseek"
            / "foldseek"
            / "bin"
            / "foldseek"
        ),
    )
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
        args.output_root or here / "confirm" / "results" / args.subset / "foldseek"
    ).resolve()

    records = load_records(jsonl_path, args.subset, structures)
    binary_hash = sha256(binary)
    version = foldseek_version(binary)
    print(
        f"method=foldseek subset={args.subset} records={len(records)} "
        f"workers={args.workers} output={output_root} "
        f"foldseek_version={version} binary_sha256={binary_hash}",
        flush=True,
    )

    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                run_one,
                record,
                binary,
                binary_hash,
                version,
                output_root,
                args.overwrite,
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
