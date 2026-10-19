#!/usr/bin/env python3
"""Run DeepAlign on unique SABmark TwiFP or SupFP structure pairs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


INPUT_RE = re.compile(
    r"^(?P<which>[12])(?:st|nd) input protein:\s*(?P<name>\S+)\s+length=\s*(?P<length>\d+)",
    re.MULTILINE,
)
NORM_RE = re.compile(r"Normalization length .* is\s+(\d+)")
DCUTOFF_RE = re.compile(r"Distance cutoff value is\s*([0-9.eE+-]+)")
SCORE_HEADER = "# BLOSUM CLESUM DeepScore SeqID LALI RMSD(A) TMscore MAXSUB GDT_TS GDT_HA uGDT"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_output(text: str):
    inputs = {int(m.group("which")): m for m in INPUT_RE.finditer(text)}
    norm = NORM_RE.search(text)
    cutoff = DCUTOFF_RE.search(text)
    lines = text.splitlines()
    try:
        header_index = next(i for i, line in enumerate(lines) if line.strip() == SCORE_HEADER)
        score_line = next(
            line.strip() for line in lines[header_index + 1:]
            if line.strip() and not line.startswith("#")
        )
    except (StopIteration, ValueError) as exc:
        raise ValueError("missing DeepAlign score record") from exc
    fields = score_line.split()
    if set(inputs) != {1, 2} or norm is None or cutoff is None or len(fields) != 11:
        raise ValueError("incomplete DeepAlign output")
    return {
        "name_a": inputs[1].group("name"),
        "name_b": inputs[2].group("name"),
        "length_a": int(inputs[1].group("length")),
        "length_b": int(inputs[2].group("length")),
        "normalization_length": int(norm.group(1)),
        "distance_cutoff": float(cutoff.group(1)),
        "blosum": float(fields[0]),
        "clesum": float(fields[1]),
        "deep_score": float(fields[2]),
        "seqid": float(fields[3]),
        "aligned_length": int(fields[4]),
        "rmsd": float(fields[5]),
        "tmscore": float(fields[6]),
        "maxsub": float(fields[7]),
        "gdt_ts": float(fields[8]),
        "gdt_ha": float(fields[9]),
        "ugdt": float(fields[10]),
    }


def metadata_path(output_root: Path, job_id: str) -> Path:
    return output_root / "jobs" / job_id[:6] / job_id / "result.json"


def write_json_atomic(path: Path, value) -> None:
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def reusable(result, destination: Path, binary_hash: str) -> bool:
    return (
        result.get("status") == "ok"
        and result.get("binary_sha256") == binary_hash
        and (destination / "stdout.txt").is_file()
        and bool(result.get("metrics"))
    )


def run_one(job, binary: Path, binary_hash: str, output_root: Path,
            timeout: int, overwrite: bool):
    job_id = job["job_id"]
    metadata = metadata_path(output_root, job_id)
    destination = metadata.parent
    destination.mkdir(parents=True, exist_ok=True)
    if metadata.is_file() and not overwrite:
        previous = json.loads(metadata.read_text())
        if reusable(previous, destination, binary_hash):
            previous["cached"] = True
            return previous

    result = {
        "job_id": job_id,
        "sid_a": job["sid_a"], "sid_b": job["sid_b"],
        "pdb_a": job["pdb_a"], "pdb_b": job["pdb_b"],
        "benchmark_occurrences": int(job["benchmark_occurrences"]),
        "binary": str(binary), "binary_sha256": binary_hash,
        "status": "error", "error": "", "cached": False,
    }
    for field in ("pdb_a", "pdb_b"):
        if not Path(job[field]).is_file():
            result["error"] = f"missing input: {job[field]}"
            write_json_atomic(metadata, result)
            return result

    command = [str(binary), job["pdb_a"], job["pdb_b"]]
    result["command"] = command
    started = time.perf_counter()
    try:
        with tempfile.TemporaryDirectory(prefix=f"deepalign_{job_id}_") as scratch:
            completed = subprocess.run(
                command, cwd=scratch, capture_output=True, text=True, timeout=timeout
            )
        result["seconds"] = round(time.perf_counter() - started, 6)
        result["returncode"] = completed.returncode
        (destination / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
        (destination / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
        if completed.returncode != 0:
            result["error"] = f"returncode={completed.returncode}"
        else:
            result["metrics"] = parse_output(completed.stdout)
            result["status"] = "ok"
    except subprocess.TimeoutExpired as exc:
        result["seconds"] = round(time.perf_counter() - started, 6)
        result["status"] = "timeout"
        result["error"] = str(exc)
    except Exception as exc:
        result["seconds"] = round(time.perf_counter() - started, 6)
        result["error"] = f"{type(exc).__name__}: {exc}"
    write_json_atomic(metadata, result)
    return result


def select_jobs(jobs, query_index: Path, queries):
    if not queries:
        return jobs
    requested = set(queries)
    rows = read_tsv(query_index)
    known = {row["query_sid"] for row in rows}
    unknown = requested - known
    if unknown:
        raise ValueError(f"unknown true query SID(s): {', '.join(sorted(unknown))}")
    selected = {row["job_id"] for row in rows if row["query_sid"] in requested}
    return [job for job in jobs if job["job_id"] in selected]


def write_summary(output_root: Path, results) -> Path:
    fields = (
        "job_id", "sid_a", "sid_b", "benchmark_occurrences", "status", "cached",
        "seconds", "name_a", "name_b", "length_a", "length_b",
        "normalization_length", "distance_cutoff", "blosum", "clesum",
        "deep_score", "seqid", "aligned_length", "rmsd", "tmscore",
        "maxsub", "gdt_ts", "gdt_ha", "ugdt", "error", "result_json",
    )
    path = output_root / "deepalign_results.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for result in sorted(results, key=lambda row: row["job_id"]):
            metrics = result.get("metrics", {})
            row = {field: result.get(field, metrics.get(field, "")) for field in fields}
            row["result_json"] = str(metadata_path(output_root, result["job_id"]))
            writer.writerow(row)
    return path


def main() -> int:
    here = Path(__file__).resolve().parent
    repo = here.parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("subset", choices=("twifp", "supfp"))
    parser.add_argument("--jobs", type=Path)
    parser.add_argument("--query-index", type=Path)
    parser.add_argument("--query", action="append", default=[], help="True query SID; repeatable")
    parser.add_argument("--binary", type=Path, default=repo / "baselines/DeepAlign/DeepAlign")
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.workers < 1 or args.timeout < 1:
        parser.error("--workers and --timeout must be positive")

    subset_root = repo / "dataset/SABmark" / args.subset
    jobs_path = args.jobs or subset_root / f"manifests/{args.subset}_ealign_jobs.tsv"
    query_index = (
        args.query_index
        or subset_root / f"manifests/{args.subset}_query_targets.tsv"
    )
    output_root = (args.output_root or subset_root / "results/deepalign").resolve()
    binary = args.binary.resolve(strict=True)
    jobs = read_tsv(jobs_path.resolve(strict=True))
    jobs = select_jobs(jobs, query_index.resolve(strict=True), args.query)
    binary_hash = sha256(binary)
    print(
        f"subset={args.subset} jobs={len(jobs)} workers={args.workers} output={output_root} "
        f"binary_sha256={binary_hash}", flush=True
    )

    selected_results = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                run_one, job, binary, binary_hash, output_root,
                args.timeout, args.overwrite,
            ): job
            for job in jobs
        }
        for index, future in enumerate(as_completed(futures), 1):
            result = future.result()
            selected_results.append(result)
            if index % 100 == 0 or index == len(jobs):
                ok = sum(row["status"] == "ok" for row in selected_results)
                cached = sum(bool(row.get("cached")) for row in selected_results)
                print(f"progress={index}/{len(jobs)} ok={ok} cached={cached}", flush=True)

    all_results = [
        json.loads(path.read_text())
        for path in (output_root / "jobs").glob("*/*/result.json")
    ]
    summary = write_summary(output_root, all_results)
    failed = [row for row in selected_results if row["status"] != "ok"]
    elapsed = sum(float(row.get("seconds", 0)) for row in selected_results)
    print(f"summary={summary}")
    print(
        f"selected_ok={len(selected_results) - len(failed)} "
        f"selected_failed={len(failed)} saved_jobs={len(all_results)} "
        f"summed_pair_seconds={elapsed:.3f}"
    )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
