#!/usr/bin/env python3
"""Score pairwise structure aligners against curated manual alignments.

MALISAM, MALIDUP and HOMSTRAD truth: only uppercase-uppercase columns.
Dali input is run_homstrad_dali_direct.py output, not released *.dali.ali.
Errors are separate; the primary all-pair macro gives them an empty prediction.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter
from pathlib import Path

from run_pairwise_alignment import (BENCHMARKS, REPO, Case,
                                    ca_sequence, load_cases, sequence, sha256)
from split_HOMSTRAD import AA3

PAIRWISE_ROOT = REPO / "test_results/pairwise_runs"
DALI_ROOTS = {
    name: PAIRWISE_ROOT / "dali" / name
    for name in BENCHMARKS
}
DEFAULT_METHODS = (
    "ealign",
    "deepalign",
    "tmalign",
    "dali",
    "gtalign_cuda",
    "gtalign_patched_cpu",
    "gtalign_patched_cpu_t20",
    "foldseek",
)
METHODS = DEFAULT_METHODS
GTALIGN_METHODS = ("gtalign_cuda", "gtalign_patched_cpu")


def is_gtalign_method(method: str) -> bool:
    return method in GTALIGN_METHODS or re.fullmatch(r"gtalign_patched_cpu_t[2-9][0-9]*", method) is not None


def is_supported_method(method: str) -> bool:
    return method in METHODS or is_gtalign_method(method)


def residue_pairs(a: str, b: str, uppercase_only: bool = False,
                  offset_a: int = 0, offset_b: int = 0) -> set[tuple[int, int]]:
    if len(a) != len(b):
        raise ValueError("Alignment rows have different lengths")
    i, j = offset_a, offset_b
    pairs = set()
    for left, right in zip(a, b):
        if left.isalpha() and right.isalpha() and (not uppercase_only or left.isupper() and right.isupper()):
            pairs.add((i, j))
        i += left.isalpha()
        j += right.isalpha()
    return pairs


def fasta_rows(path: Path) -> tuple[str, str]:
    rows, current = [], []
    for line in path.read_text(errors="replace").splitlines():
        line = line.strip()
        if line.startswith(">"):
            if current:
                rows.append("".join(current))
                current = []
        elif line:
            current.append(line)
    if current:
        rows.append("".join(current))
    if len(rows) != 2:
        raise ValueError(f"Expected two FASTA rows: {path}")
    return rows[0], rows[1]


def tmalign_rows(path: Path) -> tuple[str, str]:
    lines = path.read_text(errors="replace").splitlines()
    marker = next((i for i, line in enumerate(lines)
                   if line.startswith('(":" denotes residue pairs')), None)
    if marker is None:
        raise ValueError(f"Missing TM-align alignment block: {path}")
    block = [line.strip() for line in lines[marker + 1:] if line.strip()]
    if len(block) < 3:
        raise ValueError(f"Incomplete TM-align alignment block: {path}")
    return block[0], block[2]


def checked_local_prediction(rows: tuple[str, str], struct: tuple[str, str],
                             starts: tuple[int, int], ends: tuple[int, int],
                             source: Path) -> set[tuple[int, int]]:
    """Validate 1-based inclusive local ranges, then map aligned columns globally."""
    parsed = tuple(sequence(row) for row in rows)
    expected = tuple(seq[start - 1:end] for seq, start, end in zip(struct, starts, ends))
    if min(*starts) < 1 or any(end < start for start, end in zip(starts, ends)):
        raise ValueError(f"Invalid local alignment range in {source}: {starts}..{ends}")
    if parsed != expected:
        raise ValueError(
            f"Output/PDB sequence-slice mismatch: {source}; "
            f"starts={starts} ends={ends} parsed={tuple(map(len, parsed))} "
            f"expected={tuple(map(len, expected))}"
        )
    return residue_pairs(*rows, offset_a=starts[0] - 1, offset_b=starts[1] - 1)


def gtalign_prediction(path: Path, struct: tuple[str, str]) -> set[tuple[int, int]]:
    payload = json.loads(path.read_text())
    search = payload.get("gtalign_search")
    if not isinstance(search, dict):
        raise ValueError(f"Missing gtalign_search object: {path}")
    hits = search.get("search_results", [])
    if len(hits) != 1 or not isinstance(hits[0].get("hit_record"), dict):
        raise ValueError(f"Expected exactly one GTalign hit: {path}")
    alignment = hits[0]["hit_record"].get("alignment")
    if not isinstance(alignment, dict):
        raise ValueError(f"Missing GTalign alignment object: {path}")
    rows = (str(alignment["query_aln"]), str(alignment["refrn_aln"]))
    starts = (int(alignment["query_from"]), int(alignment["refn_from"]))
    ends = (int(alignment["query_to"]), int(alignment["refn_to"]))
    return checked_local_prediction(rows, struct, starts, ends, path)


def pdb_ca_coordinates(path: Path) -> tuple[str, list[tuple[float, float, float]]]:
    letters, coordinates = [], []
    previous = None
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or line[12:16].strip() != "CA":
            continue
        aa = AA3.get(line[17:20].strip().upper())
        if aa is None:
            continue
        key = (line[21], line[22:26], line[26])
        if key == previous:
            continue
        previous = key
        letters.append(aa)
        coordinates.append(tuple(float(line[start:start + 8]) for start in (30, 38, 46)))
    return "".join(letters), coordinates


def foldseek_ca_map(sequence_text: str, ca_text: str, pdb_path: Path,
                    source: Path) -> list[int]:
    values = [float(value) for value in ca_text.rstrip(",").split(",") if value]
    if len(values) != 3 * len(sequence_text):
        raise ValueError(f"Foldseek sequence/CA count mismatch: {source}")
    foldseek_coordinates = list(zip(values[0::3], values[1::3], values[2::3]))
    pdb_sequence, pdb_coordinates = pdb_ca_coordinates(pdb_path)
    cell_size = 0.02
    coordinate_index: dict[tuple[str, int, int, int], list[int]] = {}
    for pdb_index, (letter, xyz) in enumerate(zip(pdb_sequence, pdb_coordinates)):
        key = (letter, *(math.floor(value / cell_size) for value in xyz))
        coordinate_index.setdefault(key, []).append(pdb_index)
    mapping = []
    for position, (letter, xyz) in enumerate(zip(sequence_text.upper(), foldseek_coordinates)):
        cell = tuple(math.floor(value / cell_size) for value in xyz)
        candidates = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    key = (letter, cell[0] + dx, cell[1] + dy, cell[2] + dz)
                    for pdb_index in coordinate_index.get(key, []):
                        pdb_xyz = pdb_coordinates[pdb_index]
                        squared = sum((left - right) ** 2 for left, right in zip(xyz, pdb_xyz))
                        if squared <= 0.0004:
                            candidates.append((squared, pdb_index))
        if not candidates:
            raise ValueError(
                f"Cannot map Foldseek residue {position + 1}/{len(sequence_text)} "
                f"to {pdb_path}: {source}"
            )
        candidates.sort()
        if len(candidates) > 1 and abs(candidates[0][0] - candidates[1][0]) < 1e-12:
            raise ValueError(f"Ambiguous Foldseek CA coordinate at residue {position + 1}: {source}")
        mapping.append(candidates[0][1])
    return mapping


def mapped_residue_pairs(rows: tuple[str, str], starts: tuple[int, int],
                          ends: tuple[int, int], maps: tuple[list[int], list[int]],
                          source: Path) -> set[tuple[int, int]]:
    if len(rows[0]) != len(rows[1]):
        raise ValueError(f"Foldseek alignment rows have different lengths: {source}")
    positions = [starts[0] - 1, starts[1] - 1]
    pairs = set()
    for left, right in zip(*rows):
        mapped = [None, None]
        for side, letter in enumerate((left, right)):
            if letter.isalpha():
                if not 0 <= positions[side] < len(maps[side]):
                    raise ValueError(f"Foldseek alignment index is out of range: {source}")
                mapped[side] = maps[side][positions[side]]
                positions[side] += 1
        if mapped[0] is not None and mapped[1] is not None:
            pairs.add((mapped[0], mapped[1]))
    if tuple(positions) != ends:
        raise ValueError(f"Foldseek alignment range mismatch: {source}; consumed={tuple(positions)} ends={ends}")
    return pairs


def foldseek_prediction(path: Path, pdb_paths: tuple[Path, Path]) -> set[tuple[int, int]]:
    lines = [line for line in path.read_text(errors="replace").splitlines() if line.strip()]
    if len(lines) != 1:
        raise ValueError(f"Expected exactly one Foldseek row: {path}")
    fields = (
        "query", "target", "qstart", "qend", "tstart", "tend", "qlen", "tlen",
        "qseq", "tseq", "qca", "tca", "qaln", "taln", "bits", "evalue", "prob", "alnlen",
    )
    values = lines[0].split("\t")
    if len(values) != len(fields):
        raise ValueError(f"Expected {len(fields)} Foldseek fields, found {len(values)}: {path}")
    row = dict(zip(fields, values))
    lengths = (int(row["qlen"]), int(row["tlen"]))
    full_sequences = (row["qseq"], row["tseq"])
    if lengths != tuple(map(len, full_sequences)):
        raise ValueError(f"Foldseek qlen/tlen mismatch full sequences: {path}")
    rows = (row["qaln"], row["taln"])
    starts = (int(row["qstart"]), int(row["tstart"]))
    ends = (int(row["qend"]), int(row["tend"]))
    expected = tuple(text[start - 1:end].upper()
                     for text, start, end in zip(full_sequences, starts, ends))
    if tuple(sequence(aligned) for aligned in rows) != expected:
        raise ValueError(f"Foldseek aligned/full-sequence mismatch: {path}")
    maps = (
        foldseek_ca_map(row["qseq"], row["qca"], pdb_paths[0], path),
        foldseek_ca_map(row["tseq"], row["tca"], pdb_paths[1], path),
    )
    return mapped_residue_pairs(rows, starts, ends, maps, path)


def binary_prediction(case: Case, method: str, root: Path) -> tuple[str, set[tuple[int, int]]]:
    folder = root / method / case.benchmark / case.case_id
    result = json.loads((folder / "result.json").read_text())
    # Result bundles may be relocated after a run.  The ordered SHA-256 values
    # below are the stable identity check; absolute paths are only provenance.
    for name, path in (("pdb_a", case.pdb_a), ("pdb_b", case.pdb_b)):
        if result["sha256"][name] != sha256(path):
            raise ValueError(f"Input changed since run: {path}")
    if result["status"] != "ok":
        return result["status"], set()
    struct = (ca_sequence(case.pdb_a), ca_sequence(case.pdb_b))
    if method == "tmalign":
        rows = tmalign_rows(folder / "stdout.txt")
        parsed = tuple(sequence(row) for row in rows)
        if parsed != struct:
            raise ValueError(f"Output/PDB sequence mismatch: {folder}; {tuple(map(len, parsed))} != {tuple(map(len, struct))}")
        prediction = residue_pairs(*rows)
    elif method in ("ealign", "deepalign"):
        rows = fasta_rows(folder / f"{case.case_id}.fasta")
        parsed = tuple(sequence(row) for row in rows)
        if parsed != struct:
            raise ValueError(f"Output/PDB sequence mismatch: {folder}; {tuple(map(len, parsed))} != {tuple(map(len, struct))}")
        prediction = residue_pairs(*rows)
    elif is_gtalign_method(method):
        prediction = gtalign_prediction(folder / "gtalign.json", struct)
    elif method == "foldseek":
        prediction = foldseek_prediction(folder / "foldseek.tsv", (case.pdb_a, case.pdb_b))
    else:
        raise ValueError(f"Unsupported binary method: {method}")
    return "ok", prediction


def pdb_dali_records(path: Path) -> dict[str, list[tuple[int, tuple[float, float, float]]]]:
    """Index PDB CA records by Dali's five-column residue ID and coordinates."""
    records: dict[str, list[tuple[int, tuple[float, float, float]]]] = {}
    previous = None
    index = -1
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or line[12:16].strip() != "CA":
            continue
        if line[17:20].strip().upper() not in AA3:
            continue
        key = (line[21], line[22:26], line[26])
        if key == previous:
            continue
        previous = key
        index += 1
        number = line[22:26] if line[26] == " " else line[23:26] + line[26]
        coordinates = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
        records.setdefault(line[21] + number, []).append((index, coordinates))
    return records


def dat_map(dat_path: Path, pdb_path: Path) -> dict[int, int]:
    """Map Dali's 1-based imported positions to 0-based PDB CA positions."""
    lines = dat_path.read_text(errors="replace").splitlines()
    header = re.match(r"^>>>>\s+\S+\s+(\d+)\b", lines[0]) if lines else None
    if header is None:
        raise ValueError(f"Malformed Dali DAT: {dat_path}")
    expected = int(header[1])
    records = pdb_dali_records(pdb_path)
    mapped, used = {}, set()
    for line in lines:
        if not line.startswith("-acc "):
            continue
        position = int(line[5:9])
        key = line[11:16]
        xyz = tuple(float(line[start:start + 12]) for start in (22, 34, 46))
        choices = [(sum((x - y) ** 2 for x, y in zip(coords, xyz)), index)
                   for index, coords in records.get(key, []) if index not in used
                   and all(abs(x - y) <= 0.11 for x, y in zip(coords, xyz))]
        if not choices:
            raise ValueError(f"Cannot map Dali residue {key!r} in {dat_path} to {pdb_path}")
        index = min(choices)[1]
        mapped[position] = index
        used.add(index)
    if set(mapped) != set(range(1, expected + 1)):
        raise ValueError(f"Incomplete Dali import mapping: {dat_path}")
    return mapped


def dali_prediction(case: Case, root: Path) -> tuple[str, set[tuple[int, int]]]:
    folder = root / "pairs" / case.case_id
    result = json.loads((folder / "result.json").read_text())
    expected = (case.pdb_a.resolve(), case.pdb_b.resolve())
    expected_hashes = tuple(sha256(path) for path in expected)
    recorded_hashes = (
        result["pdb_sha256"]["query"],
        result["pdb_sha256"]["target"],
    )
    if recorded_hashes == expected_hashes:
        source = expected
        swap = False
    elif recorded_hashes == expected_hashes[::-1]:
        source = expected[::-1]
        swap = True
    else:
        raise ValueError(f"Dali input hashes do not match benchmark: {folder}")
    if result["status"] == "no_hit":
        return "no_hit", set()
    if result["status"] != "hit":
        return result["status"], set()
    dat_a = list((folder / "DAT").glob("q001*.dat"))
    dat_b = list((folder / "DAT").glob("t001*.dat"))
    if len(dat_a) != 1 or len(dat_b) != 1:
        raise ValueError(f"Expected one DAT file per Dali input: {folder}")
    map_a, map_b = dat_map(dat_a[0], source[0]), dat_map(dat_b[0], source[1])
    prediction = set()
    with (folder / "alignment_pairs.tsv").open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            i, j = int(row["query_dali_index"]), int(row["target_dali_index"])
            pair = (map_a[i], map_b[j])
            prediction.add(pair[::-1] if swap else pair)
    if not prediction:
        raise ValueError(f"Dali hit has no equivalence pairs: {folder}")
    return "hit", prediction


def score(truth: set[tuple[int, int]], prediction: set[tuple[int, int]]) -> dict:
    tp = len(truth & prediction)
    precision = tp / len(prediction) if prediction else 0.0
    recall = tp / len(truth) if truth else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {"truth_pairs": len(truth), "predicted_pairs": len(prediction), "tp": tp,
            "precision": precision, "recall": recall, "f1": f1}


def evaluate(case: Case, method: str, binary_root: Path, dali_roots: dict[str, Path]) -> dict:
    output = (dali_roots[case.benchmark] / "pairs" / case.case_id if method == "dali"
              else binary_root / method / case.benchmark / case.case_id)
    row = {"method": method, "benchmark": case.benchmark, "id": case.case_id,
           "pdb_a": str(case.pdb_a), "pdb_b": str(case.pdb_b),
           "manual": str(case.manual), "output": str(output), "status": "error", "error": ""}
    truth = residue_pairs(case.row_a, case.row_b, uppercase_only=True)
    if not truth:
        raise ValueError(f"No uppercase-uppercase reference pairs: {case.manual}")
    try:
        if method == "dali":
            status, prediction = dali_prediction(case, dali_roots[case.benchmark])
        else:
            status, prediction = binary_prediction(case, method, binary_root)
        row["status"] = status
        if status in ("ok", "hit", "no_hit"):
            row.update(score(truth, prediction))
        else:
            saved = json.loads((output / "result.json").read_text())
            row["error"] = saved.get("error") or f"run status={status}"
    except (OSError, ValueError, KeyError, IndexError, json.JSONDecodeError) as exc:
        row["status"] = "evaluation_error"
        row["error"] = str(exc)
    return row


def summarize(rows: list[dict]) -> dict:
    scored = [row for row in rows if row["status"] in ("ok", "hit", "no_hit")]
    all_values = lambda key: sum(float(row.get(key, 0.0)) for row in rows) / len(rows)
    scored_values = lambda key: sum(float(row[key]) for row in scored) / len(scored) if scored else None
    return {"pairs": len(rows), "statuses": dict(Counter(row["status"] for row in rows)),
            "coverage": len(scored) / len(rows),
            "macro_precision": all_values("precision"), "macro_recall": all_values("recall"),
            "macro_f1": all_values("f1"),
            "bad_cases_f1_lt_0.5": sum(float(row.get("f1", 0.0)) < 0.5 for row in rows),
            "successful_only_macro_precision": scored_values("precision"),
            "successful_only_macro_recall": scored_values("recall"),
            "successful_only_macro_f1": scored_values("f1")}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmarks", nargs="+", choices=BENCHMARKS, default=list(BENCHMARKS))
    parser.add_argument(
        "--methods", nargs="+", default=list(DEFAULT_METHODS),
        help="Methods to score; multithreaded CPU GTalign uses gtalign_patched_cpu_tN",
    )
    parser.add_argument("--binary-root", type=Path, default=PAIRWISE_ROOT)
    parser.add_argument("--dali-homstrad", type=Path, default=DALI_ROOTS["homstrad"])
    parser.add_argument("--dali-malisam", type=Path, default=DALI_ROOTS["malisam"])
    parser.add_argument("--dali-malidup", type=Path, default=DALI_ROOTS["malidup"])
    parser.add_argument("--split-root", type=Path, default=REPO / "test_results/HOMSTRAD_split_structures")
    parser.add_argument(
        "--output-root",
        type=Path,
        default=REPO / "test_results/01_pairwise_alignment/evaluation",
    )
    parser.add_argument("--limit", type=int, default=0, help="First N pairs per benchmark")
    parser.add_argument("--allow-partial", action="store_true", help="Permit missing or invalid result files in a smoke test")
    args = parser.parse_args()
    unsupported = [method for method in args.methods if not is_supported_method(method)]
    if unsupported:
        parser.error(f"unsupported --methods: {', '.join(unsupported)}")
    args.binary_root = args.binary_root.resolve()
    args.output_root = args.output_root.resolve()
    cases = load_cases(args.benchmarks, args.split_root)
    if args.limit:
        cases = [case for benchmark in args.benchmarks
                 for case in [c for c in cases if c.benchmark == benchmark][:args.limit]]
    dali_roots = {name: getattr(args, f"dali_{name}").resolve() for name in BENCHMARKS}
    rows = [evaluate(case, method, args.binary_root, dali_roots)
            for method in args.methods for case in cases]
    rows.sort(key=lambda row: (row["benchmark"], row["method"], row["id"]))
    args.output_root.mkdir(parents=True, exist_ok=True)
    detail = args.output_root / "per_pair.csv"
    columns = ("benchmark", "method", "id", "status", "truth_pairs", "predicted_pairs",
               "tp", "precision", "recall", "f1", "pdb_a", "pdb_b", "manual", "output", "error")
    with detail.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in columns} for row in rows)
    groups = {(benchmark, method): summarize([row for row in rows
              if row["benchmark"] == benchmark and row["method"] == method])
              for benchmark in args.benchmarks for method in args.methods}
    summary = args.output_root / "summary.csv"
    fields = ("benchmark", "method", "pairs", "ok", "hit", "no_hit", "error", "timeout", "evaluation_error",
              "coverage", "macro_precision", "macro_recall", "macro_f1", "bad_cases_f1_lt_0.5",
              "successful_only_macro_precision", "successful_only_macro_recall", "successful_only_macro_f1")
    with summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for (benchmark, method), value in groups.items():
            writer.writerow({"benchmark": benchmark, "method": method, **
                             {key: value["statuses"].get(key, 0) for key in ("ok", "hit", "no_hit", "error", "timeout", "evaluation_error")},
                             **{key: value[key] for key in fields if key in value}})
    (args.output_root / "summary.json").write_text(json.dumps(
        {f"{benchmark}/{method}": value for (benchmark, method), value in groups.items()}, indent=2) + "\n")
    errors = [row for row in rows if row["status"] not in ("ok", "hit", "no_hit")]
    print(f"pairs={len(cases)} methods={len(args.methods)} errors={len(errors)}")
    print(f"per_pair={detail}\nsummary={summary}")
    if any(row["status"] == "evaluation_error" for row in rows) and not args.allow_partial:
        raise SystemExit("Missing/invalid result files; inspect per_pair.csv or use --allow-partial for a smoke test")


if __name__ == "__main__":
    main()
