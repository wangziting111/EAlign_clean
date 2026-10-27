#!/usr/bin/env python3
"""Run one pairwise aligner and save its raw outputs without scoring.

Examples:
  python3 experiments/01_pairwise_alignment/run_pairwise_alignment.py --method ealign --benchmarks malisam --limit 1
  python3 experiments/01_pairwise_alignment/run_pairwise_alignment.py --method gtalign_cuda --workers 1 --limit 1
  python3 experiments/01_pairwise_alignment/run_pairwise_alignment.py --method gtalign_patched_cpu --workers 1 --cpu-threads 20 --limit 1
  python3 experiments/01_pairwise_alignment/run_pairwise_alignment.py --method foldseek --workers 1 --limit 1
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

REPO = Path(__file__).resolve().parents[2]
DATASET = REPO / "dataset"
SPLIT = REPO / "test_results/HOMSTRAD_split_structures"
DEFAULT_OUTPUT = REPO / "test_results/pairwise_runs"
DEFAULT_BINARIES = {
    "ealign": REPO / "binary/ealign",
    "deepalign": REPO / "binary/baselines/DeepAlign/DeepAlign",
    "tmalign": REPO / "binary/baselines/TMalign/TMalign",
    "gtalign_cuda": (
        REPO / "binary/baselines/gtalign_alpha/Linux_installer_GPU/bin/gtalign"
    ),
    "gtalign_patched_cpu": (
        REPO / "binary/baselines/gtalign_alpha/patched_cpu/bin/gtalign"
    ),
    "foldseek": REPO / "binary/baselines/foldseek/bin/foldseek",
}
GTALIGN_METHODS = ("gtalign_cuda", "gtalign_patched_cpu")
FOLDSEEK_FIELDS = (
    "query", "target", "qstart", "qend", "tstart", "tend", "qlen", "tlen",
    "qseq", "tseq", "qca", "tca", "qaln", "taln", "bits", "evalue", "prob", "alnlen",
)
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


def run_configuration(method: str, cpu_threads: int, benchmark: str) -> dict:
    if method in GTALIGN_METHODS:
        configuration = {
            "hetatm": True,
            "dev_min_length": 3,
            "speed": 0,
            "pre_score": 0,
            "report_score_threshold": 0,
            "add_search_by_ss": True,
            "dev_mem_mb": 4096,
            "nhits": 1,
            "nalns": 1,
            "output_format": "json",
        }
        if method == "gtalign_patched_cpu":
            configuration.update({"cpu_threads": cpu_threads, "cpu_threads_reading": 1})
        return configuration
    if method == "foldseek":
        configuration = {
            "alignment_type": 2,
            "exhaustive_search": True,
            "evalue": "inf",
            "max_seqs": 1,
            "threads": 1,
            "format_output": ",".join(FOLDSEEK_FIELDS),
        }
        if benchmark == "malidup":
            configuration["pdb_cleanup"] = True
        return configuration
    return {}


def output_method(method: str, cpu_threads: int) -> str:
    """Keep the 1-thread baseline name and isolate every multithreaded CPU run."""
    if method == "gtalign_patched_cpu" and cpu_threads > 1:
        return f"gtalign_patched_cpu_t{cpu_threads}"
    return method


def expected_alignment(destination: Path, method: str, case_id: str) -> Path:
    if method == "tmalign":
        return destination / "stdout.txt"
    if method in ("ealign", "deepalign"):
        return destination / f"{case_id}.fasta"
    if method in GTALIGN_METHODS:
        return destination / "gtalign.json"
    if method == "foldseek":
        return destination / "foldseek.tsv"
    raise ValueError(f"Unsupported method: {method}")


def clean_malidup_pdb(source: Path, destination: Path) -> None:
    """Remove MALIDUP's legacy coordinate-record fields after column 66."""
    lines = source.read_text(errors="replace").splitlines()
    normalized = [line[:66] if line.startswith(("ATOM  ", "HETATM")) else line
                  for line in lines]
    destination.write_text("\n".join(normalized) + "\n")


def build_command(case: Case, method: str, binary: Path, scratch: Path,
                  cpu_threads: int) -> tuple[list[str], Path | None]:
    if method == "tmalign":
        return [str(binary), str(case.pdb_a.resolve()), str(case.pdb_b.resolve())], None
    if method in ("ealign", "deepalign"):
        return ([str(binary), str(case.pdb_a.resolve()), str(case.pdb_b.resolve()),
                 "-p", "1", "-o", case.case_id], None)
    if method in GTALIGN_METHODS:
        output = scratch / "gtalign_output"
        output.mkdir()
        command = [
            str(binary), f"--qrs={case.pdb_a.resolve()}", f"--rfs={case.pdb_b.resolve()}",
            "-o", str(output), "--hetatm", "--dev-min-length=3", "--speed=0",
            "--pre-score=0", "-s", "0", "--add-search-by-ss", "--dev-mem=4096",
            "--nhits=1", "--nalns=1", "--outfmt=1",
        ]
        if method == "gtalign_patched_cpu":
            command += [f"--cpu-threads={cpu_threads}", "--cpu-threads-reading=1"]
        return command, output
    if method == "foldseek":
        query = case.pdb_a.resolve()
        target = case.pdb_b.resolve()
        if case.benchmark == "malidup":
            query = scratch / "query.pdb"
            target = scratch / "target.pdb"
            clean_malidup_pdb(case.pdb_a, query)
            clean_malidup_pdb(case.pdb_b, target)
        output = scratch / "foldseek.tsv"
        command = [
            str(binary), "easy-search", str(query), str(target),
            str(output), str(scratch / "foldseek_tmp"), "--alignment-type", "2",
            "--exhaustive-search", "1", "-e", "inf", "--max-seqs", "1",
            "--threads", "1", "--format-output", ",".join(FOLDSEEK_FIELDS),
        ]
        return command, output
    raise ValueError(f"Unsupported method: {method}")


def run_one(case: Case, method: str, binary: Path, binary_hash: str,
            method_label: str, output_root: Path, timeout: int,
            overwrite: bool, cpu_threads: int) -> dict:
    destination = output_root / method_label / case.benchmark / case.case_id
    metadata = destination / "result.json"
    input_hashes = {"pdb_a": sha256(case.pdb_a), "pdb_b": sha256(case.pdb_b), "binary": binary_hash}
    configuration = run_configuration(method, cpu_threads, case.benchmark)
    if metadata.is_file() and not overwrite:
        saved = json.loads(metadata.read_text())
        if saved.get("sha256") != input_hashes:
            raise ValueError(f"Inputs/binary changed for {metadata}; use --overwrite or a new output root")
        if saved.get("timing_scope") != TIMING_SCOPE:
            raise ValueError(f"Cached timing scope is obsolete for {metadata}; use --overwrite or a new output root")
        if method in (*GTALIGN_METHODS, "foldseek") and saved.get("run_configuration") != configuration:
            raise ValueError(f"Run configuration changed for {metadata}; use --overwrite or a new output root")
        if saved.get("status") == "ok":
            alignment = expected_alignment(destination, method, case.case_id)
            if not alignment.is_file():
                raise FileNotFoundError(f"Cached result is missing its alignment: {alignment}")
        return saved
    destination.mkdir(parents=True, exist_ok=True)
    result = {
        "benchmark": case.benchmark, "id": case.case_id, "method": method_label,
        "pdb_a": str(case.pdb_a.resolve()), "pdb_b": str(case.pdb_b.resolve()),
        "manual": str(case.manual.resolve()), "binary": str(binary),
        "sha256": input_hashes, "timing_scope": TIMING_SCOPE,
        "run_configuration": configuration,
        "status": "error", "error": "",
    }
    with tempfile.TemporaryDirectory(prefix=f"{method}_{case.case_id}_") as scratch:
        scratch_path = Path(scratch)
        command, raw_output = build_command(case, method, binary, scratch_path, cpu_threads)
        result["command"] = command
        started = time.perf_counter()
        try:
            done = subprocess.run(command, cwd=scratch, capture_output=True, text=True, timeout=timeout)
            result["seconds"] = round(time.perf_counter() - started, 3)
            result["returncode"] = done.returncode
            (destination / "stdout.txt").write_text(done.stdout)
            (destination / "stderr.txt").write_text(done.stderr)
            if method in ("ealign", "deepalign"):
                for suffix in (".fasta", ".score", ".local", ".pdb"):
                    source = scratch_path / f"{case.case_id}{suffix}"
                    if source.exists():
                        shutil.copy2(source, destination / source.name)
                if done.returncode == 0 and (scratch_path / f"{case.case_id}.fasta").is_file():
                    result["status"] = "ok"
            elif method == "tmalign" and done.returncode == 0 and done.stdout.strip():
                result["status"] = "ok"
            elif method in GTALIGN_METHODS and done.returncode == 0:
                json_files = sorted(raw_output.glob("*.json")) if raw_output else []
                if len(json_files) == 1:
                    shutil.copy2(json_files[0], destination / "gtalign.json")
                    result["status"] = "ok"
                else:
                    result["error"] = f"expected one GTalign JSON output, found {len(json_files)}"
            elif method == "foldseek" and done.returncode == 0 and raw_output and raw_output.is_file():
                rows = [line for line in raw_output.read_text(errors="replace").splitlines() if line.strip()]
                if len(rows) == 1 and len(rows[0].split("\t")) == len(FOLDSEEK_FIELDS):
                    shutil.copy2(raw_output, destination / "foldseek.tsv")
                    result["status"] = "ok"
                else:
                    result["error"] = f"expected one {len(FOLDSEEK_FIELDS)}-field Foldseek row, found {len(rows)} rows"
            if result["status"] != "ok":
                if not result["error"]:
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
    parser.add_argument("--cpu-threads", type=int, default=1,
                        help="Threads inside each gtalign_patched_cpu process (default: 1)")
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--limit", type=int, default=0, help="First N pairs per benchmark")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    if args.cpu_threads < 1:
        parser.error("--cpu-threads must be positive")
    args.output_root = args.output_root.resolve()
    cases = load_cases(args.benchmarks, args.split_root)
    if args.limit:
        cases = [case for benchmark in args.benchmarks
                 for case in [c for c in cases if c.benchmark == benchmark][:args.limit]]
    binary = (args.binary or DEFAULT_BINARIES[args.method]).resolve(strict=True)
    binary_hash = sha256(binary)
    method_label = output_method(args.method, args.cpu_threads)
    print(f"method={method_label} pairs={len(cases)} workers={args.workers}", flush=True)
    if args.method == "gtalign_cuda" and args.workers > 1:
        print("warning: use --workers 1 unless multiple GTalign processes are intentionally sharing the GPU", flush=True)
    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_one, case, args.method, binary, binary_hash,
                               method_label, args.output_root, args.timeout,
                               args.overwrite, args.cpu_threads): case
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
