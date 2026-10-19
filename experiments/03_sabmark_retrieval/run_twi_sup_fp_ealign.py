#!/usr/bin/env python3
"""Run EAlign on unique SABmark TwiFP or SupFP pairs, with resume support."""

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
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


ARTIFACT_SUFFIXES = (".fasta", ".score", ".local", ".pdb")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_score(path: Path):
    for line in path.read_text(errors="replace").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 21 or fields[4] != "->" or fields[8] != "->":
            raise ValueError(f"unexpected score record: {line}")
        return {
            "name1": fields[0], "name2": fields[1],
            "len1": int(fields[2]), "len2": int(fields[3]),
            "blosum": float(fields[5]), "score_3di": float(fields[6]),
            "combined_score": float(fields[7]), "lali": int(fields[9]),
            "rmsd": float(fields[10]), "tmscore": float(fields[11]),
            "maxsub": float(fields[13]), "gdt_ts": float(fields[14]),
            "gdt_ha": float(fields[15]), "seqid": float(fields[17]),
            "normalization_length": int(fields[18]),
            "dcutoff": float(fields[19]), "ugdt": float(fields[20]),
        }
    raise ValueError(f"no score record in {path}")


def metadata_path(output_root: Path, job_id: str) -> Path:
    return output_root / "jobs" / job_id[:6] / job_id / "result.json"


def result_is_reusable(result, destination: Path, binary_hash: str) -> bool:
    if result.get("status") != "ok" or result.get("binary_sha256") != binary_hash:
        return False
    return all((destination / f"{result['job_id']}{suffix}").is_file()
               for suffix in ARTIFACT_SUFFIXES)


def write_json_atomic(path: Path, value) -> None:
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def run_one(job, binary: Path, binary_hash: str, output_root: Path,
            timeout: int, overwrite: bool):
    job_id = job["job_id"]
    metadata = metadata_path(output_root, job_id)
    destination = metadata.parent
    destination.mkdir(parents=True, exist_ok=True)
    if metadata.is_file() and not overwrite:
        previous = json.loads(metadata.read_text())
        if result_is_reusable(previous, destination, binary_hash):
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

    with tempfile.TemporaryDirectory(prefix=f"ealign_{job_id}_") as scratch_name:
        scratch = Path(scratch_name)
        command = [str(binary), job["pdb_a"], job["pdb_b"], "-p", "1", "-o", job_id]
        result["command"] = command
        started = time.perf_counter()
        try:
            done = subprocess.run(
                command, cwd=scratch, capture_output=True, text=True, timeout=timeout
            )
            result["seconds"] = round(time.perf_counter() - started, 6)
            result["returncode"] = done.returncode
            (destination / "stdout.txt").write_text(done.stdout, encoding="utf-8")
            (destination / "stderr.txt").write_text(done.stderr, encoding="utf-8")
            produced = [scratch / f"{job_id}{suffix}" for suffix in ARTIFACT_SUFFIXES]
            if done.returncode != 0:
                result["error"] = f"returncode={done.returncode}"
            elif not all(path.is_file() for path in produced):
                result["error"] = "missing one or more expected EAlign output files"
            else:
                for source in produced:
                    shutil.move(str(source), destination / source.name)
                result["metrics"] = parse_score(destination / f"{job_id}.score")
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
    job_ids = {row["job_id"] for row in rows if row["query_sid"] in requested}
    return [job for job in jobs if job["job_id"] in job_ids]


def write_summary(output_root: Path, results) -> Path:
    fields = (
        "job_id", "sid_a", "sid_b", "benchmark_occurrences", "status", "cached",
        "seconds", "combined_score", "tmscore", "lali", "rmsd", "maxsub",
        "gdt_ts", "gdt_ha", "seqid", "ugdt", "error", "result_json",
    )
    path = output_root / "ealign_results.tsv"
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
    parser.add_argument("--binary", type=Path, default=repo / "binary/ealign")
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
    output_root = (args.output_root or subset_root / "results/ealign").resolve()
    binary = args.binary.resolve(strict=True)
    jobs = read_tsv(jobs_path.resolve(strict=True))
    jobs = select_jobs(jobs, query_index.resolve(strict=True), args.query)
    binary_hash = sha256(binary)
    print(
        f"subset={args.subset} jobs={len(jobs)} workers={args.workers} output={output_root} "
        f"binary_sha256={binary_hash}", flush=True
    )

    results = []
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
            results.append(result)
            if index % 100 == 0 or index == len(jobs):
                ok = sum(row["status"] == "ok" for row in results)
                cached = sum(bool(row.get("cached")) for row in results)
                print(f"progress={index}/{len(jobs)} ok={ok} cached={cached}", flush=True)

    all_results = [
        json.loads(path.read_text())
        for path in (output_root / "jobs").glob("*/*/result.json")
    ]
    summary = write_summary(output_root, all_results)
    failed = [row for row in results if row["status"] != "ok"]
    elapsed = sum(float(row.get("seconds", 0)) for row in results)
    print(f"summary={summary}")
    print(
        f"selected_ok={len(results) - len(failed)} selected_failed={len(failed)} "
        f"saved_jobs={len(all_results)} summed_pair_seconds={elapsed:.3f}"
    )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
