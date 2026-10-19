#!/usr/bin/env python3
"""Run one pairwise aligner and save its raw outputs without scoring.

Examples:
  python test_config/run_pairwise_alignment.py --method ealign --benchmarks malisam --limit 1
  python test_config/run_pairwise_alignment.py --method tmalign --workers 4
Results: OUTPUT/METHOD/BENCHMARK/CASE/{result.json, alignment files, logs}.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from split_HOMSTRAD import HOMSTRAD, ca_sequence, parse_ali, sequence

REPO = Path(__file__).resolve().parents[1]
DATASET = REPO / "dataset"
SPLIT = REPO / "test_results/HOMSTRAD_split_structures"
DEFAULT_OUTPUT = REPO / "test_results/pairwise_runs"
DEFAULT_BINARIES = {
    "ealign": REPO / "binary/ealign",
    "deepalign": REPO / "baselines/DeepAlign/DeepAlign",
    "tmalign": REPO / "baselines/TMalign/TMalign",
}
BENCHMARKS = ("malisam", "malidup", "homstrad")
EXPECTED = {"malisam": 130, "malidup": 241, "homstrad": 1032}
TIMING_SCOPE = "binary_subprocess_wall_time"


@dataclass(frozen=True)
class Case:
    benchmark: str
    case_id: str
    pdb_a: Path
    pdb_b: Path
    manual: Path
    row_a: str
    row_b: str


def read_manual(path: Path) -> tuple[str, str]:
    rows = [line.strip() for line in path.read_text(errors="replace").splitlines()
            if line.strip() and not line.startswith(">")]
    if len(rows) != 2 or len(rows[0]) != len(rows[1]):
        raise ValueError(f"Expected two equal-length manual rows: {path}")
    return rows[0], rows[1]


def load_cases(benchmarks: list[str], split_root: Path = SPLIT) -> list[Case]:
    cases = []
    for benchmark in benchmarks:
        current = []
        if benchmark in ("malisam", "malidup"):
            root = DATASET / ("Malisam" if benchmark == "malisam" else "Malidup")
            pair_list = DATASET / ("pdblistsam.csv" if benchmark == "malisam" else "pdblistdup.csv")
            for record in csv.DictReader(pair_list.open(newline="")):
                case_id = record["dir"]
                manuals = list((root / case_id).glob("*.manual.ali"))
                if len(manuals) != 1:
                    raise ValueError(f"Expected one manual alignment for {case_id}: {manuals}")
                row_a, row_b = read_manual(manuals[0])
                paths = [root / record[name] for name in ("file1", "file2")]
                if any(not p.is_file() for p in paths):
                    raise FileNotFoundError(paths)
                seqs = [ca_sequence(p) for p in paths]
                references = [sequence(row_a), sequence(row_b)]
                if seqs == references:
                    pass
                elif seqs == references[::-1]:
                    paths.reverse()  # Always run binaries in manual-row order.
                else:
                    raise ValueError(f"PDB/manual sequence mismatch: {case_id}")
                current.append(Case(benchmark, case_id, paths[0], paths[1], manuals[0], row_a, row_b))
        else:
            for manual in sorted(HOMSTRAD.glob("*/*.ali")):
                family = manual.parent.name
                order, rows = parse_ali(manual)
                first, last = order[0], order[-1]
                pdb_a = split_root / family / f"{first}.pdb"
                pdb_b = split_root / family / f"{last}.pdb"
                if not pdb_a.is_file() or not pdb_b.is_file():
                    raise FileNotFoundError(f"Run split_HOMSTRAD.py first: {family}")
                for pdb, row in ((pdb_a, rows[first]), (pdb_b, rows[last])):
                    struct, ref = ca_sequence(pdb), sequence(row)
                    if not ref.startswith(struct) or len(ref) - len(struct) > 1:
                        raise ValueError(f"Split/manual sequence mismatch: {pdb}")
                current.append(Case(benchmark, family, pdb_a, pdb_b, manual, rows[first], rows[last]))
        if len(current) != EXPECTED[benchmark] or len({c.case_id for c in current}) != len(current):
            raise ValueError(f"Expected {EXPECTED[benchmark]} unique {benchmark} pairs, got {len(current)}")
        cases.extend(sorted(current, key=lambda c: c.case_id))
    return cases


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_one(case: Case, method: str, binary: Path, binary_hash: str,
            output_root: Path, timeout: int, overwrite: bool) -> dict:
    destination = output_root / method / case.benchmark / case.case_id
    metadata = destination / "result.json"
    input_hashes = {"pdb_a": sha256(case.pdb_a), "pdb_b": sha256(case.pdb_b), "binary": binary_hash}
    if metadata.is_file() and not overwrite:
        saved = json.loads(metadata.read_text())
        if saved.get("sha256") != input_hashes:
            raise ValueError(f"Inputs/binary changed for {metadata}; use --overwrite or a new output root")
        if saved.get("timing_scope") != TIMING_SCOPE:
            raise ValueError(f"Cached timing scope is obsolete for {metadata}; use --overwrite or a new output root")
        if saved.get("status") == "ok":
            alignment = destination / ("stdout.txt" if method == "tmalign" else f"{case.case_id}.fasta")
            if not alignment.is_file():
                raise FileNotFoundError(f"Cached result is missing its alignment: {alignment}")
        return saved
    destination.mkdir(parents=True, exist_ok=True)
    result = {
        "benchmark": case.benchmark, "id": case.case_id, "method": method,
        "pdb_a": str(case.pdb_a.resolve()), "pdb_b": str(case.pdb_b.resolve()),
        "manual": str(case.manual.resolve()), "binary": str(binary),
        "sha256": input_hashes, "timing_scope": TIMING_SCOPE,
        "status": "error", "error": "",
    }
    with tempfile.TemporaryDirectory(prefix=f"{method}_{case.case_id}_") as scratch:
        command = [str(binary), str(case.pdb_a.resolve()), str(case.pdb_b.resolve())]
        if method != "tmalign":
            command += ["-p", "1", "-o", case.case_id]
        result["command"] = command
        started = time.perf_counter()
        try:
            done = subprocess.run(command, cwd=scratch, capture_output=True, text=True, timeout=timeout)
            result["seconds"] = round(time.perf_counter() - started, 3)
            result["returncode"] = done.returncode
            (destination / "stdout.txt").write_text(done.stdout)
            (destination / "stderr.txt").write_text(done.stderr)
            if method != "tmalign":
                for suffix in (".fasta", ".score", ".local", ".pdb"):
                    source = Path(scratch) / f"{case.case_id}{suffix}"
                    if source.exists():
                        shutil.copy2(source, destination / source.name)
                if done.returncode == 0 and (Path(scratch) / f"{case.case_id}.fasta").is_file():
                    result["status"] = "ok"
            elif done.returncode == 0 and done.stdout.strip():
                result["status"] = "ok"
            if result["status"] != "ok":
                result["error"] = f"returncode={done.returncode}; missing expected alignment output"
        except subprocess.TimeoutExpired as exc:
            result["seconds"] = round(time.perf_counter() - started, 3)
            result["status"] = "timeout"
            result["error"] = str(exc)
    metadata.write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmarks", nargs="+", choices=BENCHMARKS, default=list(BENCHMARKS))
    parser.add_argument("--method", required=True, choices=tuple(DEFAULT_BINARIES))
    parser.add_argument("--binary", type=Path, help="Override the selected method's default binary")
    parser.add_argument("--split-root", type=Path, default=SPLIT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--limit", type=int, default=0, help="First N pairs per benchmark")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    args.output_root = args.output_root.resolve()
    cases = load_cases(args.benchmarks, args.split_root)
    if args.limit:
        cases = [case for benchmark in args.benchmarks
                 for case in [c for c in cases if c.benchmark == benchmark][:args.limit]]
    binary = (args.binary or DEFAULT_BINARIES[args.method]).resolve(strict=True)
    binary_hash = sha256(binary)
    print(f"method={args.method} pairs={len(cases)} workers={args.workers}", flush=True)
    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_one, case, args.method, binary, binary_hash,
                               args.output_root, args.timeout, args.overwrite): case
                   for case in cases}
        for index, future in enumerate(as_completed(futures), 1):
            rows.append(future.result())
            if index % 25 == 0 or index == len(cases):
                print(f"{index}/{len(cases)} ok={sum(r['status'] == 'ok' for r in rows)}", flush=True)
    # Summarize all saved runs, including earlier --limit invocations.
    all_results = [json.loads(path.read_text()) for path in args.output_root.glob("*/*/*/result.json")]
    summary = args.output_root / "run_summary.csv"
    summary.parent.mkdir(parents=True, exist_ok=True)
    with summary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("method", "benchmark", "id", "status", "seconds", "timing_scope", "error"),
        )
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in writer.fieldnames}
                         for row in sorted(all_results, key=lambda row: (row["method"], row["benchmark"], row["id"])))
    print(f"summary={summary}")


if __name__ == "__main__":
    main()
