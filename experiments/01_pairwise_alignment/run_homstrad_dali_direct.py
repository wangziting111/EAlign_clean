#!/usr/bin/env python3
"""Run DaliLite v5 pairwise on HOMSTRAD, MALISAM, or MALIDUP.

MALISAM/MALIDUP pairs come from pdblistsam.csv/pdblistdup.csv. If a PDB
lacks HEADER, only the per-pair temporary Dali input gets one. Outputs
are under OUTPUT/pairs/CASE/: the original Dali report, DAT files,
logs, Dali-index pairs, and result.json. This script does not score
against the released reference alignments.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import tempfile
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path


EALIGN = Path(__file__).resolve().parents[2]
HOMSTRAD = EALIGN / "data" / "HOMSTRAD"
SPLIT = EALIGN / "test_results/HOMSTRAD_split_structures"
DATASET = EALIGN / "dataset"
DALI = EALIGN / "baselines/DaliLite.v5/bin/dali.pl"
EXPECTED_PAIRS = {"homstrad": 1032, "malisam": 130, "malidup": 241}
MALI_INPUTS = {
    "malisam": (DATASET / "pdblistsam.csv", DATASET / "Malisam"),
    "malidup": (DATASET / "pdblistdup.csv", DATASET / "Malidup"),
}
OUTFMT = "summary,alignments,equivalences,transrot"
FIELDS = (
    "benchmark", "family", "query", "target", "ca_len_a", "ca_len_b", "imported_len_a",
    "imported_len_b", "status", "z_score", "rmsd", "lali",
    "equivalence_segments", "seconds", "error",
)
HIT = re.compile(
    r"^\s*1:\s+(\S+)\s+(-?\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+"
    r"(\d+)\s+\d+\s+\d+", re.MULTILINE
)
SEGMENT = re.compile(
    r"^\s*1:\s+\S+\s+\S+\s+(\d+)\s*-\s*(\d+)\s*<=>\s*"
    r"(\d+)\s*-\s*(\d+)", re.MULTILINE
)


@dataclass(frozen=True)
class Pair:
    benchmark: str
    family: str
    query: str
    target: str
    pdb_a: Path
    pdb_b: Path
    ca_len_a: int
    ca_len_b: int


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ca_count(path: Path) -> int:
    with path.open(errors="replace") as handle:
        count = 0
        previous_id = None
        chains = set()
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")) or line[12:16].strip() != "CA":
                continue
            residue_id = (line[21], line[22:26], line[26])
            chains.add(line[21])
            if residue_id != previous_id:
                count += 1
                previous_id = residue_id
    if not count or len(chains) != 1:
        raise ValueError(f"Expected one nonempty CA chain: {path}")
    return count


def load_homstrad_pairs(homstrad: Path, split: Path) -> list[Pair]:
    alignments = sorted(homstrad.glob("*/*.ali"))
    if len(alignments) != 1032:
        raise ValueError(f"Expected 1032 HOMSTRAD alignments, found {len(alignments)}")
    pairs = []
    for ali in alignments:
        members = [line[4:].strip() for line in ali.read_text().splitlines()
                   if line.startswith(">P1;")]
        if len(members) < 2 or len(members) != len(set(members)):
            raise ValueError(f"Invalid member list: {ali}")
        family, query, target = ali.parent.name, members[0], members[-1]
        pdb_a = (split / family / f"{query}.pdb").resolve(strict=True)
        pdb_b = (split / family / f"{target}.pdb").resolve(strict=True)
        pairs.append(Pair("homstrad", family, query, target, pdb_a, pdb_b,
                          ca_count(pdb_a), ca_count(pdb_b)))
    if len({pair.family for pair in pairs}) != len(pairs):
        raise ValueError("More than one alignment found for a family")
    return pairs


def load_mali_pairs(benchmark: str) -> list[Pair]:
    pair_list, root = MALI_INPUTS[benchmark]
    root = root.resolve(strict=True)
    with pair_list.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not {"dir", "file1", "file2"} <= set(reader.fieldnames or []):
            raise ValueError(f"Expected dir,file1,file2 columns: {pair_list}")
        rows = list(reader)
    if len(rows) != EXPECTED_PAIRS[benchmark]:
        raise ValueError(f"Expected {EXPECTED_PAIRS[benchmark]} {benchmark} pairs, found {len(rows)}")
    pairs = []
    for row in rows:
        family = row["dir"]
        paths = []
        for field in ("file1", "file2"):
            relative = Path(row[field])
            if relative.is_absolute() or relative.parent != Path(family):
                raise ValueError(f"Unexpected {field} path for {family}: {relative}")
            path = (root / relative).resolve(strict=True)
            path.relative_to(root)  # Refuse paths escaping the dataset root.
            paths.append(path)
        pdb_a, pdb_b = paths
        pairs.append(Pair(benchmark, family, pdb_a.stem, pdb_b.stem,
                          pdb_a, pdb_b, ca_count(pdb_a), ca_count(pdb_b)))
    if len({pair.family for pair in pairs}) != len(pairs):
        raise ValueError(f"Duplicate case IDs in {pair_list}")
    return pairs


def dat_length(path: Path) -> int:
    with path.open(errors="replace") as handle:
        first = handle.readline()
    match = re.match(r"^>>>>\s+\S+\s+(\d+)\b", first)
    if not match or int(match[1]) == 0:
        raise ValueError(f"Empty or malformed Dali DAT: {path}")
    return int(match[1])


def parse_report(report: str, target_id: str) -> tuple[dict, list[tuple[int, int]]]:
    if "# No:" not in report or "# Structural equivalences" not in report:
        raise ValueError("Incomplete Dali report")
    summary, equivalences = report.split("# Structural equivalences", 1)
    hit = HIT.search(summary)
    segments = SEGMENT.findall(equivalences)
    if hit is None:
        if segments:
            raise ValueError("Structural equivalences without a summary hit")
        return {"status": "no_hit", "lali": 0, "equivalence_segments": 0}, []
    if hit[1].replace("-", "") != target_id:
        raise ValueError(f"Unexpected Dali target: {hit[1]}")
    residue_pairs = []
    for qa, qb, ta, tb in segments:
        qa, qb, ta, tb = map(int, (qa, qb, ta, tb))
        if qb < qa or tb - ta != qb - qa:
            raise ValueError(f"Invalid equivalence segment: {qa}-{qb} <=> {ta}-{tb}")
        residue_pairs.extend((qa + offset, ta + offset)
                             for offset in range(qb - qa + 1))
    lali = int(hit[4])
    if len(residue_pairs) != lali:
        raise ValueError(f"Dali lali={lali}, parsed {len(residue_pairs)} pairs")
    return {
        "status": "hit", "z_score": float(hit[2]), "rmsd": float(hit[3]),
        "lali": lali, "equivalence_segments": len(segments),
    }, residue_pairs


def run_dali(command: list[str], cwd: Path, timeout: int) -> tuple[int, str, str, float, bool]:
    start = time.perf_counter()
    process = subprocess.Popen(command, cwd=cwd, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, text=True, errors="replace",
                               start_new_session=True)
    timed_out = False
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        stdout, stderr = process.communicate()
    return process.returncode, stdout, stderr, round(time.perf_counter() - start, 3), timed_out


def save_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def prepare_pdb(source: Path, destination: Path) -> bool:
    """Use the original PDB if headed; otherwise add HEADER only in scratch."""
    with source.open("rb") as original:
        first = original.readline()
        if first.startswith(b"HEADER"):
            destination.symlink_to(source)
            return False
        with destination.open("wb") as prepared:
            prepared.write(b"HEADER   \n")
            prepared.write(first)
            shutil.copyfileobj(original, prepared)
    return True


def run_pair(pair: Pair, dali: Path, output: Path, timeout: int,
             tool_hashes: dict[str, str]) -> dict:
    pair_dir = output / "pairs" / pair.family
    pair_dir.mkdir(parents=True, exist_ok=True)
    result_path = pair_dir / "result.json"
    pdb_hashes = {"query": sha256(pair.pdb_a), "target": sha256(pair.pdb_b)}
    if result_path.exists():
        result = json.loads(result_path.read_text())
        if (result.get("benchmark", "homstrad") != pair.benchmark
                or result.get("pdb_sha256") != pdb_hashes
                or result.get("tool_sha256") != tool_hashes
                or result.get("outfmt") != OUTFMT):
            raise ValueError(f"Inputs or Dali build changed: {pair_dir}; use a new output directory")
        result.setdefault("benchmark", pair.benchmark)
        return result

    command = [
        str(dali), "--pdbfile1", "q.pdb", "--pdbid1", "q001",
        "--pdbfile2", "t.pdb", "--pdbid2", "t001",
        "--dat1", "DAT", "--dat2", "DAT",
        "--outfmt", OUTFMT, "--title", pair.family,
    ]
    result = {
        "benchmark": pair.benchmark, "family": pair.family,
        "query": pair.query, "target": pair.target,
        "pdb_a": str(pair.pdb_a), "pdb_b": str(pair.pdb_b),
        "ca_len_a": pair.ca_len_a, "ca_len_b": pair.ca_len_b,
        "pdb_sha256": pdb_hashes, "tool_sha256": tool_hashes,
        "command": command, "outfmt": OUTFMT, "status": "error", "error": "",
    }
    with tempfile.TemporaryDirectory(prefix=f"dali_{pair.benchmark}_") as scratch:
        tmp = Path(scratch)
        (tmp / "DAT").mkdir()
        result["header_added_a"] = prepare_pdb(pair.pdb_a, tmp / "q.pdb")
        result["header_added_b"] = prepare_pdb(pair.pdb_b, tmp / "t.pdb")
        code, stdout, stderr, seconds, timed_out = run_dali(command, tmp, timeout)
        (pair_dir / "stdout.log").write_text(stdout)
        (pair_dir / "stderr.log").write_text(stderr)
        result.update(returncode=code, seconds=seconds)
        for dat in (tmp / "DAT").glob("*.dat"):
            destination = pair_dir / "DAT"
            destination.mkdir(exist_ok=True)
            shutil.copy2(dat, destination / dat.name)
        reports = list(tmp.glob("q001*.txt"))
        for report in reports:
            shutil.copy2(report, pair_dir / report.name)
        try:
            if timed_out:
                raise TimeoutError(f"Exceeded {timeout} seconds")
            if code != 0:
                raise ValueError(f"dali.pl exit code {code}")
            query_dat = list((tmp / "DAT").glob("q001*.dat"))
            target_dat = list((tmp / "DAT").glob("t001*.dat"))
            if len(query_dat) != 1 or len(target_dat) != 1 or len(reports) != 1:
                raise ValueError("Expected one query DAT, target DAT, and Dali report")
            result["imported_len_a"] = dat_length(query_dat[0])
            result["imported_len_b"] = dat_length(target_dat[0])
            result["report"] = reports[0].name
            parsed, residue_pairs = parse_report(
                reports[0].read_text(errors="replace"), target_dat[0].stem
            )
            result.update(parsed)
            with (pair_dir / "alignment_pairs.tsv").open("w") as handle:
                handle.write("query_dali_index\ttarget_dali_index\n")
                handle.writelines(f"{a}\t{b}\n" for a, b in residue_pairs)
        except (OSError, ValueError, TimeoutError) as exc:
            result.update(status="timeout" if timed_out else "error", error=str(exc))
    save_json(result_path, result)
    return result


def save_summary(path: Path, rows: dict[str, dict]) -> None:
    temporary = path.with_suffix(".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows[family] for family in sorted(rows))
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", choices=tuple(EXPECTED_PAIRS), default="homstrad")
    parser.add_argument("--homstrad-root", type=Path, default=HOMSTRAD)
    parser.add_argument("--split-root", type=Path, default=SPLIT)
    parser.add_argument("--dali-bin", type=Path, default=DALI)
    parser.add_argument("--output-root", type=Path, help="Default: test_results/<benchmark>_dali")
    parser.add_argument("--family", action="append", help="Run only this family/case ID; repeatable")
    parser.add_argument("--limit", type=int, default=0, help="Run the first N selected pairs")
    parser.add_argument("--timeout", type=int, default=600, help="Seconds per pair")
    parser.add_argument("--jobs", type=int, default=1, help="Independent Dali pairs to run at once")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    if args.limit < 0 or args.timeout <= 0 or args.jobs <= 0:
        parser.error("--limit must be nonnegative; --timeout and --jobs positive")
    pairs = (load_homstrad_pairs(args.homstrad_root, args.split_root)
             if args.benchmark == "homstrad" else load_mali_pairs(args.benchmark))
    if args.family:
        selected = set(args.family)
        missing = selected - {pair.family for pair in pairs}
        if missing:
            parser.error(f"Unknown families: {sorted(missing)}")
        pairs = [pair for pair in pairs if pair.family in selected]
    if args.limit:
        pairs = pairs[:args.limit]
    print(f"Preflight: {EXPECTED_PAIRS[args.benchmark]} {args.benchmark} pairs; "
          f"selected {len(pairs)}", flush=True)
    short = [pair.family for pair in pairs if min(pair.ca_len_a, pair.ca_len_b) <= 29]
    if short:
        print(f"Warning: {len(short)} selected pairs have a chain of at most 29 CA "
              f"residues; Dali may skip them: {', '.join(short)}", flush=True)
    if args.check_only:
        return 0

    dali = args.dali_bin.resolve(strict=True)
    tools = {"dali.pl": dali, **{
        name: dali.parent / name for name in ("mpidali.pm", "mkdssp", "serialcompare")
    }}
    tool_hashes = {name: sha256(path) for name, path in tools.items()}
    output = (args.output_root or EALIGN / "test_results" / f"{args.benchmark}_dali").resolve()
    output.mkdir(parents=True, exist_ok=True)
    rows = {}
    for saved in sorted((output / "pairs").glob("*/result.json")):
        result = json.loads(saved.read_text())
        if (result.get("benchmark", "homstrad") != args.benchmark
                or result.get("tool_sha256") != tool_hashes
                or result.get("outfmt") != OUTFMT):
            raise ValueError(f"Different benchmark or Dali build already recorded: {saved}")
        result.setdefault("benchmark", args.benchmark)
        rows[result["family"]] = result
    selected_rows = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(run_pair, pair, dali, output, args.timeout, tool_hashes): pair
            for pair in pairs
        }
        for index, future in enumerate(as_completed(futures), 1):
            pair = futures[future]
            result = future.result()
            selected_rows.append(result)
            rows[pair.family] = result
            save_summary(output / "summary.csv", rows)
            print(f"[{index}/{len(pairs)}] {pair.family}: {result['status']} "
                  f"Z={result.get('z_score', '')} Lali={result.get('lali', '')} "
                  f"t={result.get('seconds', '')}s", flush=True)
    counts = Counter(row["status"] for row in selected_rows)
    print(f"Finished: {dict(counts)}; {output / 'summary.csv'}", flush=True)
    return 1 if counts["error"] or counts["timeout"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
