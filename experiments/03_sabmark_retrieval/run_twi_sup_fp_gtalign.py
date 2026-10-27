#!/usr/bin/env python3
"""Run GTalign CPU or GPU once for every TT/TF SABmark JSONL record."""

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
BASE_RUN_CONFIGURATION = {
    "hetatm": True,
    "dev_min_length": 3,
    "speed": 0,
    "pre_score": 0,
    "report_score_threshold": 0,
    "add_search_by_ss": True,
    "dev_mem_mb": 4096,
    "output_format": "json",
}
INPUT_INTERPRETATION = "gtalign_default_ter3_split0_with_hetatm"


def run_configuration(backend: str, cpu_threads: int) -> dict:
    configuration = dict(BASE_RUN_CONFIGURATION)
    # Keep the historical GPU configuration byte-for-byte compatible so a
    # paused GPU run can reuse its existing result.json files.
    if backend == "cpu":
        configuration.update(
            {
                "backend": "cpu",
                "cpu_threads": cpu_threads,
            }
        )
    return configuration


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def inspect_gtalign_binary(binary: Path) -> tuple[str, str]:
    completed = subprocess.run(
        [str(binary), "-h"],
        capture_output=True,
        text=True,
    )
    help_text = completed.stdout + "\n" + completed.stderr
    version = ""
    for line in help_text.splitlines():
        if line.lower().startswith("gtalign "):
            version = line.strip()
            break
    if not version:
        raise RuntimeError("could not determine GTalign version from -h output")

    if "CPU/multiprocessing compilation without support for GPUs" in help_text:
        backend = "cpu"
    elif "compiled with GPU support" in help_text:
        backend = "gpu"
    else:
        raise RuntimeError("could not determine whether GTalign is a CPU or GPU build")
    return version, backend


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


def parse_gtalign_json(path: Path) -> tuple[str, dict]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    search = payload.get("gtalign_search")
    if not isinstance(search, dict):
        raise ValueError("missing gtalign_search object")

    hits = search.get("search_results", [])
    if len(hits) != 1:
        raise ValueError(f"expected exactly one GTalign hit, found {len(hits)}")
    hit = hits[0].get("hit_record")
    if not isinstance(hit, dict):
        raise ValueError("missing hit_record object")
    alignment = hit.get("alignment")
    if not isinstance(alignment, dict):
        raise ValueError("missing alignment object")

    query_length = int(hit["query_length"])
    reference_length = int(hit["reference_length"])
    tm_query = float(alignment["tmscore_query"])
    tm_reference = float(alignment["tmscore_refn"])

    if query_length <= reference_length:
        tm_minlen, tm_maxlen = tm_query, tm_reference
    else:
        tm_minlen, tm_maxlen = tm_reference, tm_query
    tm_hmean = (
        0.0
        if tm_query + tm_reference == 0.0
        else 2.0 * tm_query * tm_reference / (tm_query + tm_reference)
    )

    metrics = {
        "query_description": search.get("query", {}).get("description", ""),
        "reference_description": hit.get("reference_description", ""),
        "query_length": query_length,
        "reference_length": reference_length,
        "tm_score_query": tm_query,
        "tm_score_reference": tm_reference,
        "tm_score_greater": max(tm_query, tm_reference),
        "tm_score_minlen": tm_minlen,
        "tm_score_maxlen": tm_maxlen,
        "tm_score_hmean": tm_hmean,
        "rmsd": float(alignment["rmsd"]),
        "n_matched": int(alignment["n_matched"]),
        "n_aligned": int(alignment["n_aligned"]),
        "n_identities": int(alignment["n_identities"]),
        "n_gaps": int(alignment["n_gaps"]),
        "query_from": int(alignment["query_from"]),
        "query_to": int(alignment["query_to"]),
        "reference_from": int(alignment["refn_from"]),
        "reference_to": int(alignment["refn_to"]),
    }
    return str(search.get("version", "")), metrics


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
    configuration: dict,
) -> bool:
    return (
        result.get("status") == "ok"
        and result.get("binary_sha256") == binary_hash
        and result.get("run_configuration") == configuration
        and result.get("input_interpretation") == INPUT_INTERPRETATION
        and result.get("pdb1") == record["pdb1"]
        and result.get("pdb2") == record["pdb2"]
        and (destination / "gtalign.json").is_file()
        and bool(result.get("metrics"))
    )


def run_one(
    record: dict,
    binary: Path,
    binary_hash: str,
    version_text: str,
    backend: str,
    cpu_threads: int,
    configuration: dict,
    output_root: Path,
    overwrite: bool,
) -> dict:
    destination = record_directory(output_root, record)
    destination.mkdir(parents=True, exist_ok=True)
    result_path = destination / "result.json"

    if result_path.is_file() and not overwrite:
        previous = json.loads(result_path.read_text())
        if reusable(
            previous, destination, record, binary_hash, configuration
        ):
            previous["cached"] = True
            return previous

    result = {
        **record,
        "binary": str(binary),
        "binary_sha256": binary_hash,
        "gtalign_build": version_text,
        "backend": backend,
        "cpu_threads": cpu_threads if backend == "cpu" else None,
        "run_configuration": configuration,
        "input_interpretation": INPUT_INTERPRETATION,
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
            prefix=f"gtalign_{backend}_{record['record_id']}_"
        ) as scratch_name:
            scratch = Path(scratch_name)
            gtalign_output = scratch / "output"
            gtalign_output.mkdir()
            command = [
                str(binary),
                f"--qrs={record['pdb1']}",
                f"--rfs={record['pdb2']}",
                "-o",
                str(gtalign_output),
                "--hetatm",
                "--dev-min-length=3",
                "--speed=0",
                "--pre-score=0",
                "-s",
                "0",
                "--add-search-by-ss",
                "--dev-mem=4096",
                "--outfmt=1",
            ]
            if backend == "cpu":
                command.append(f"--cpu-threads={cpu_threads}")
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
            else:
                json_files = sorted(gtalign_output.glob("*.json"))
                if len(json_files) != 1:
                    result["error"] = (
                        "expected exactly one GTalign JSON file, "
                        f"found {len(json_files)}"
                    )
                else:
                    gtalign_version_value, metrics = parse_gtalign_json(
                        json_files[0]
                    )
                    result["gtalign_version"] = gtalign_version_value
                    result["metrics"] = metrics
                    shutil.copy2(json_files[0], destination / "gtalign.json")
                    result["status"] = "ok"
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"

    write_json_atomic(result_path, result)
    return result


def write_summary(
    output_root: Path,
    results: list[dict],
    backend: str,
    cpu_threads: int,
) -> Path:
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
        "backend",
        "cpu_threads",
        "alignment_seconds",
        "gtalign_version",
        "query_length",
        "reference_length",
        "tm_score_query",
        "tm_score_reference",
        "tm_score_greater",
        "tm_score_minlen",
        "tm_score_maxlen",
        "tm_score_hmean",
        "rmsd",
        "n_matched",
        "n_aligned",
        "n_identities",
        "n_gaps",
        "query_from",
        "query_to",
        "reference_from",
        "reference_to",
        "input_interpretation",
        "error",
        "result_json",
    )
    path = output_root / "gtalign_results.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for result in sorted(results, key=lambda row: row["source_line"]):
            metrics = result.get("metrics", {})
            row = {
                field: result.get(field, metrics.get(field, ""))
                for field in fields
            }
            # Historical cached GPU records predate the backend field.
            row["backend"] = backend
            row["cpu_threads"] = cpu_threads if backend == "cpu" else ""
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
        "--backend",
        choices=("gpu", "cpu"),
        default="gpu",
        help="GTalign build to use (default: gpu)",
    )
    parser.add_argument(
        "--binary",
        type=Path,
        help="override the binary selected by --backend",
    )
    parser.add_argument(
        "--cpu-threads",
        type=int,
        default=1,
        help="threads inside each CPU GTalign process (default: 1)",
    )
    parser.add_argument("--output-root", type=Path)
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="concurrent GTalign processes (default: 1)",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    if args.workers < 1:
        parser.error("--workers must be positive")
    if not 1 <= args.cpu_threads <= 1024:
        parser.error("--cpu-threads must be in [1, 1024]")

    jsonl_path = (
        args.input_jsonl or here / args.subset / JSONL_NAMES[args.subset]
    ).resolve(strict=True)
    structures = args.structures.resolve(strict=True)
    default_binary = (
        structure_alignment_root
        / "gtalign_alpha"
        / ("Linux_installer_GPU" if args.backend == "gpu" else "Linux_installer_mp")
        / "bin"
        / "gtalign"
    )
    binary = (args.binary or default_binary).resolve(strict=True)
    output_root = (
        args.output_root
        or here / "confirm" / "results" / args.subset / f"gtalign_{args.backend}"
    ).resolve()

    records = load_records(jsonl_path, args.subset, structures)
    binary_hash = sha256(binary)
    version_text, detected_backend = inspect_gtalign_binary(binary)
    if detected_backend != args.backend:
        parser.error(
            f"--backend={args.backend} but {binary} is a {detected_backend} build"
        )
    configuration = run_configuration(args.backend, args.cpu_threads)
    print(
        f"method=gtalign_{args.backend} subset={args.subset} "
        f"records={len(records)} workers={args.workers} "
        f"cpu_threads={args.cpu_threads if args.backend == 'cpu' else 'n/a'} "
        f"output={output_root} "
        f"build={version_text!r} binary_sha256={binary_hash}",
        flush=True,
    )
    if args.backend == "gpu" and args.workers > 1:
        print(
            "warning: multiple workers share GTalign's default GPU selection; "
            "use with care",
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
                version_text,
                args.backend,
                args.cpu_threads,
                configuration,
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

    summary = write_summary(
        output_root, results, args.backend, args.cpu_threads
    )
    failed = [row for row in results if row["status"] != "ok"]
    elapsed = sum(
        float(row.get("alignment_seconds", 0.0))
        for row in results
        if row["status"] == "ok"
    )
    print(f"summary={summary}")
    print(
        f"ok={len(results) - len(failed)} failed={len(failed)} "
        f"summed_pair_seconds={elapsed:.3f}"
    )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
