#!/usr/bin/env python3
"""Score EAlign/DeepAlign/TM-align/Dali against curated manual alignments.

MALISAM, MALIDUP and HOMSTRAD truth: only uppercase-uppercase columns.
Dali input is run_homstrad_dali_direct.py output, not released *.dali.ali.
Errors are separate; the primary all-pair macro gives them an empty prediction.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path

from run_pairwise_alignment import (BENCHMARKS, DEFAULT_OUTPUT, REPO, Case,
                                    ca_sequence, load_cases, sequence, sha256)
from split_HOMSTRAD import AA3

DALI_ROOTS = {name: REPO / f"test_results/{name}_dali" for name in BENCHMARKS}
METHODS = ("ealign", "deepalign", "tmalign", "dali")


def residue_pairs(a: str, b: str, uppercase_only: bool = False) -> set[tuple[int, int]]:
    if len(a) != len(b):
        raise ValueError("Alignment rows have different lengths")
    i = j = 0
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


def binary_prediction(case: Case, method: str, root: Path) -> tuple[str, set[tuple[int, int]]]:
    folder = root / method / case.benchmark / case.case_id
    result = json.loads((folder / "result.json").read_text())
    if (Path(result["pdb_a"]).resolve(), Path(result["pdb_b"]).resolve()) != (case.pdb_a.resolve(), case.pdb_b.resolve()):
        raise ValueError(f"Input order changed: {folder}")
    for name, path in (("pdb_a", case.pdb_a), ("pdb_b", case.pdb_b)):
        if result["sha256"][name] != sha256(path):
            raise ValueError(f"Input changed since run: {path}")
    if result["status"] != "ok":
        return result["status"], set()
    rows = tmalign_rows(folder / "stdout.txt") if method == "tmalign" else fasta_rows(folder / f"{case.case_id}.fasta")
    struct = (ca_sequence(case.pdb_a), ca_sequence(case.pdb_b))
    parsed = tuple(sequence(row) for row in rows)
    if parsed != struct:
        raise ValueError(f"Output/PDB sequence mismatch: {folder}; {tuple(map(len, parsed))} != {tuple(map(len, struct))}")
    return "ok", residue_pairs(*rows)


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
    source = (Path(result["pdb_a"]).resolve(), Path(result["pdb_b"]).resolve())
    expected = (case.pdb_a.resolve(), case.pdb_b.resolve())
    if source == expected:
        swap = False
    elif source == expected[::-1]:
        swap = True
    else:
        raise ValueError(f"Dali input paths do not match benchmark: {folder}")
    for label, path in (("query", source[0]), ("target", source[1])):
        if result["pdb_sha256"][label] != sha256(path):
            raise ValueError(f"Dali input changed since run: {path}")
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
    parser.add_argument("--methods", nargs="+", choices=METHODS, default=list(METHODS))
    parser.add_argument("--binary-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--dali-homstrad", type=Path, default=DALI_ROOTS["homstrad"])
    parser.add_argument("--dali-malisam", type=Path, default=DALI_ROOTS["malisam"])
    parser.add_argument("--dali-malidup", type=Path, default=DALI_ROOTS["malidup"])
    parser.add_argument("--split-root", type=Path, default=REPO / "test_results/HOMSTRAD_split_structures")
    parser.add_argument("--output-root", type=Path, default=REPO / "test_results/pairwise_evaluation")
    parser.add_argument("--limit", type=int, default=0, help="First N pairs per benchmark")
    parser.add_argument("--allow-partial", action="store_true", help="Permit missing or invalid result files in a smoke test")
    args = parser.parse_args()
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
