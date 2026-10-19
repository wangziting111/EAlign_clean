#!/usr/bin/env python3
"""Run and evaluate an EAlign-compatible binary on the curated benchmarks.

The script keeps outputs in a separate result directory and maps
FASTA alignments back to benchmark sequence coordinates before
computing residue-pair precision/recall/F1.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from statistics import mean, median

import numpy as np


REPO = Path(__file__).resolve().parents[2]
MALI_ROOT = REPO / "dataset"
MALI_PDB_LISTS = {
    "malisam": MALI_ROOT / "pdblistsam.csv",
    "malidup": MALI_ROOT / "pdblistdup.csv",
}
HOMSTRAD_DATA = REPO / "data" / "HOMSTRAD"
DEFAULT_DEEPALIGN = REPO / "binary" / "ealign"
DEFAULT_OUTPUT_ROOT = REPO / "test_results" / "full_alignment"
HOMSTRAD_CHAIN_Dir = REPO / "test_results" / "HOMSTRAD_split_structures"
AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "ASX": "B", "GLX": "Z", "SEC": "U", "PYL": "O", "MSE": "M",
}

# Audited exceptions: the PDB chain is missing only the final .ali residue.
HOMSTRAD_MISSING_C_TERMINUS = {
    ("cat", "1m85a"),
    ("hpr", "1ka5a"),
}


def read_alignment(path: Path) -> tuple[str, str]:
    lines = [
        line.strip()
        for line in path.read_text(errors="replace").splitlines()
        if line.strip() and not line.startswith(">")
    ]
    if len(lines) != 2:
        raise ValueError(f"Expected two alignment rows in {path}, got {len(lines)}")
    if len(lines[0]) != len(lines[1]):
        raise ValueError(f"Alignment rows differ in length: {path}")
    return lines[0], lines[1]


def alignment_pairs(
    aligned_a: str,
    aligned_b: str,
    uppercase_only: bool = False,
) -> set[tuple[int, int]]:
    i = j = 0
    pairs: set[tuple[int, int]] = set()
    for aa, bb in zip(aligned_a, aligned_b):
        if (
            aa.isalpha()
            and bb.isalpha()
            and (not uppercase_only or (aa.isupper() and bb.isupper()))
        ):
            pairs.add((i, j))
        if aa.isalpha():
            i += 1
        if bb.isalpha():
            j += 1
    return pairs


def metrics(truth: set[tuple[int, int]], pred: set[tuple[int, int]]) -> dict:
    tp = len(truth & pred)
    fp = len(pred - truth)
    fn = len(truth - pred)
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "sensitivity": recall,
        "f1": f1,
        "matches": len(pred),
    }



def pdb_ca_sequence(path: Path) -> str:
    sequence: list[str] = []
    previous_residue_id = None
    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if line[12:16].strip() != "CA":
                continue
            residue_name = line[17:20].strip().upper()
            if residue_name not in AA3_TO_1:
                continue
            residue_id = (line[21], line[22:26], line[26])
            if residue_id == previous_residue_id:
                continue
            previous_residue_id = residue_id
            sequence.append(AA3_TO_1[residue_name])
    return "".join(sequence)

# ============================================================
# structure residue index -> benchmark reference residue index
# 暂不能处理内部 residue deletion / insertion
# 只允许 suffix mismatch
# ============================================================

def sequence_index_map(structure_seq: str, reference_seq: str) -> dict[int, int]:
    # pdb sequence 比 reference 长
    if structure_seq.startswith(reference_seq):
        mapped_length = len(reference_seq)
    # reference sequence 比 pdb 长 (ali mannual稍微多一节)
    elif reference_seq.startswith(structure_seq):
        mapped_length = len(structure_seq)
    else:
        raise ValueError(
            "Structure/reference sequences differ internally: "
            f"structure={len(structure_seq)}, reference={len(reference_seq)}"
        )
    return {index: index for index in range(mapped_length)}


def parse_homstrad_ali(path: Path) -> tuple[list[str], dict[str, str]]:
    order: list[str] = []
    aligned: dict[str, str] = {}
    current = None
    skip_description = False
    for raw in path.read_text(errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("C;"):
            continue
        if line.startswith(">P1;"):
            current = line[4:].strip()# >P1;1jn1a current = 1jn1a
            order.append(current)
            aligned[current] = ""
            skip_description = True
        elif current is not None and skip_description:
            skip_description = False
        elif current is not None:
            aligned[current] += line.replace("*", "")
    return order, aligned


def reference_sequence(aligned: str) -> str:
    return "".join(char.upper() for char in aligned if char.isalpha())


def pdb_chain_sequences(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    previous_ids: dict[str, tuple[str, str]] = {}
    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if line[12:16].strip() != "CA":
                continue
            residue_name = line[17:20].strip().upper()
            if residue_name not in AA3_TO_1:
                continue
            chain = line[21]
            residue_id = (line[22:26], line[26])
            if previous_ids.get(chain) == residue_id:
                continue
            previous_ids[chain] = residue_id
            sequences.setdefault(chain, []).append(AA3_TO_1[residue_name])
    return {chain: "".join(sequence) for chain, sequence in sequences.items()}

def require_unique_chain(
    candidates: list[str],
    sup_pdb: Path,
    member: str,
    match_type: str,
) -> str:
    if len(candidates) != 1:
        raise ValueError(
            f"Expected exactly one {match_type} chain for HOMSTRAD member {member} "
            f"in {sup_pdb}, found {len(candidates)}: {candidates}. "
            "Inspect this family manually instead of selecting heuristically."
        )
    return candidates[0]


def select_homstrad_chains(
    sup_pdb: Path,
    references: dict[str, str],
) -> tuple[dict[str, str], dict[str, int]]:
    chain_sequences = pdb_chain_sequences(sup_pdb)
    family = sup_pdb.parent.name
    chain_map: dict[str, str] = {}
    keep_lengths: dict[str, int] = {}

    for member, reference in references.items():
        candidates = [
            chain
            for chain, sequence in chain_sequences.items()
            if sequence == reference
        ]
        if candidates:
            chain = require_unique_chain(
                candidates, sup_pdb, member, "exact-sequence"
            )
            keep_length = len(reference)
        else:
            # Audited preprocessing rule: ignore a unique C-terminal extension.
            candidates = [
                chain
                for chain, sequence in chain_sequences.items()
                if len(sequence) > len(reference) and sequence.startswith(reference)
            ]
            if candidates:
                chain = require_unique_chain(
                    candidates, sup_pdb, member, "C-terminal-extension"
                )
                keep_length = len(reference)
            elif (family, member) in HOMSTRAD_MISSING_C_TERMINUS:
                candidates = [
                    chain
                    for chain, sequence in chain_sequences.items()
                    if len(reference) == len(sequence) + 1
                    and reference.startswith(sequence)
                ]
                chain = require_unique_chain(
                    candidates,
                    sup_pdb,
                    member,
                    "audited-missing-C-terminus",
                )
                keep_length = len(reference) - 1
            else:
                raise ValueError(
                    f"No exact-sequence chain for HOMSTRAD member {member} in "
                    f"{sup_pdb}; no audited terminal exception applies. "
                    "Inspect this family manually."
                )
        chain_map[member] = chain
        keep_lengths[member] = keep_length

    if len(set(chain_map.values())) != len(chain_map):
        raise ValueError(
            f"Selected HOMSTRAD members map to the same chain in {sup_pdb}: "
            f"{chain_map}. Inspect this family manually."
        )
    return chain_map, keep_lengths

def write_trimmed_chain(
    sup_pdb: Path,
    output_path: Path,
    chain: str,
    keep_residues: int,
    expected_sequence: str,
    overwrite: bool,
) -> None:
    if (
        output_path.exists()
        and not overwrite
        and pdb_ca_sequence(output_path) == expected_sequence
    ):
        return
    output: list[str] = []
    previous_residue_id = None
    residue_index = -1
    with sup_pdb.open(errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")) or line[21] != chain:
                continue
            residue_name = line[17:20].strip().upper()
            if residue_name not in AA3_TO_1:
                continue
            residue_id = (line[22:26], line[26])
            if residue_id != previous_residue_id:
                residue_index += 1
                previous_residue_id = residue_id
            if residue_index < keep_residues:
                output.append(line)
    if residue_index + 1 < keep_residues:
        raise ValueError(
            f"Only found {residue_index + 1} residues for chain {chain} in {sup_pdb}; "
            f"expected {keep_residues}"
        )
    header = (
        f"HEADER   \n"
    )
    output_path.write_text(header+"".join(output) + "END\n")

def prepare_homstrad_structures(
    sup_pdb: Path,
    references: dict[str, str],
    output_dir: Path,
) -> tuple[dict[str, Path], dict[str, str]]:
    chain_map, keep_lengths = select_homstrad_chains(sup_pdb, references)
    chain_sequences = pdb_chain_sequences(sup_pdb)
    paths = {member: output_dir / f"{member}.pdb" for member in references}
    for member, path in paths.items():
        keep_length = keep_lengths[member]
        write_trimmed_chain(
            sup_pdb,
            path,
            chain_map[member],
            keep_length,
            chain_sequences[chain_map[member]][:keep_length],
            overwrite=False,
        )
    return paths, chain_map


def parse_deepalign_fasta(path: Path) -> tuple[str, str]:
    records: list[str] = []
    current: list[str] = []
    for raw in path.read_text(errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current:
                records.append("".join(current))
                current = []
        else:
            current.append(line)
    if current:
        records.append("".join(current))
    if len(records) != 2:
        raise ValueError(f"Expected two FASTA records in {path}, got {len(records)}")
    if len(records[0]) != len(records[1]):
        raise ValueError(f"DeepAlign FASTA rows differ in length: {path}")
    return records[0], records[1]

# ============================================================
# 判断已有 cache 是否对应当前两个 PDB
# ============================================================

def cached_alignment_matches_structures(
    fasta_path: Path,
    pdb_a: Path,
    pdb_b: Path,
) -> bool:
    try:
        aligned_a, aligned_b = parse_deepalign_fasta(fasta_path)
        return (
            reference_sequence(aligned_a) == pdb_ca_sequence(pdb_a)
            and reference_sequence(aligned_b) == pdb_ca_sequence(pdb_b)
        )
    except (OSError, ValueError):
        return False

# ============================================================
# DeepAlign FASTA -> benchmark residue-pair coordinates
# ============================================================

def predicted_pairs_from_fasta(
    fasta_path: Path,
    map_a: dict[int, int],
    map_b: dict[int, int],
) -> set[tuple[int, int]]:
    aligned_a, aligned_b = parse_deepalign_fasta(fasta_path)
    i = j = 0
    pairs: set[tuple[int, int]] = set()
    for aa, bb in zip(aligned_a, aligned_b):
        if aa.isalpha() and bb.isalpha() and i in map_a and j in map_b:
            pairs.add((map_a[i], map_b[j]))
        if aa.isalpha():
            i += 1
        if bb.isalpha():
            j += 1
    return pairs


def parse_deepalign_score(path: Path) -> dict:
    data_line = ""
    for raw in path.read_text(errors="replace").splitlines():
        if raw.strip() and not raw.startswith("#"):
            data_line = raw.strip()
            break
    if not data_line:
        return {}
    blocks = [block.split() for block in data_line.split("->")]
    result = {"score_line": data_line}
    if len(blocks) >= 4:
        result.update({
            "name1": blocks[0][0],
            "name2": blocks[0][1],
            "len1": int(blocks[0][2]),
            "len2": int(blocks[0][3]),
            "blosum": float(blocks[1][0]),
            "clesum": float(blocks[1][1]),
            "deep_score": float(blocks[1][2]),
            "aligned_length": int(blocks[2][0]),
            "rmsd": float(blocks[2][1]),
            "tm_score": float(blocks[2][2]),
            "maxsub": float(blocks[3][0]),
            "gdt_ts": float(blocks[3][1]),
            "gdt_ha": float(blocks[3][2]),
        })
        if len(blocks) >= 5:
            result["seq_id_count"] = float(blocks[4][0])
        else:
            raise ValueError(f"missing block in .score file")
    return result

# real entry
def run_deepalign(
    binary: Path,
    pdb_a: Path,
    pdb_b: Path,
    base: str,
    output_dir: Path,
    overwrite: bool,
    reuse_only: bool,
  
) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    fasta = output_dir / f"{base}.fasta"
    local = output_dir / f"{base}.local"
    score = output_dir / f"{base}.score"
    pdb_out = output_dir / f"{base}.pdb"
    stdout = output_dir / f"{base}.stdout.txt"
    stderr = output_dir / f"{base}.stderr.txt"
    if (
        fasta.exists()
        and score.exists()
        and not overwrite
        and cached_alignment_matches_structures(fasta, pdb_a, pdb_b)
    ):
        return {"status": "cached", "fasta": fasta, "score": score}
    if reuse_only:
        return {
            "status": "error",
            "error": "cached DeepAlign FASTA/score missing or incompatible with current PDBs",
        }

    with tempfile.TemporaryDirectory(prefix=f"deepalign_{base}_") as tmp_name:
        tmp = Path(tmp_name)
        command = [
            str(binary), str(pdb_a), str(pdb_b), "-p", "1", "-o", base,
        ]
        completed = subprocess.run(
            command,
            cwd=tmp,
            capture_output=True,
            text=True,
        )# command: DeepAlign pdb_a pdb_b -p 1 -o case_id(-p output detaily)
        stdout.write_text(completed.stdout)
        stderr.write_text(completed.stderr)
        if completed.returncode != 0:
            return {
                "status": "error",
                "returncode": completed.returncode,
                "stderr_tail": completed.stderr[-1000:],
            }
        for suffix, target in {
            ".fasta": fasta,
            ".local": local,
            ".score": score,
            ".pdb": pdb_out,
        }.items():
            source = tmp / f"{base}{suffix}"
            if source.exists():
                shutil.copy2(source, target)
    if not fasta.exists() or not score.exists():
        return {"status": "error", "error": "missing DeepAlign FASTA or score"}
    return {"status": "ok", "fasta": fasta, "score": score}


def mali_pdb_pairs(dataset: str) -> dict[str, tuple[Path, Path]]:
    list_path = MALI_PDB_LISTS[dataset]
    pairs: dict[str, tuple[Path, Path]] = {}
    with list_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            case_id = row["dir"]
            if case_id in pairs:
                raise ValueError(f"Duplicate {case_id} in {list_path}")
            pdb_a = MALI_ROOT / dataset / row["file1"]
            pdb_b = MALI_ROOT / dataset / row["file2"]
            if not pdb_a.is_file() or not pdb_b.is_file():
                raise FileNotFoundError(f"Missing listed PDB for {case_id}: {pdb_a}, {pdb_b}")
            if pdb_a == pdb_b:
                raise ValueError(f"Identical listed PDB paths for {case_id}: {pdb_a}")
            pairs[case_id] = (pdb_a, pdb_b)
    return pairs


def order_listed_mali_pdbs(
    case_id: str,
    listed_a: Path,
    listed_b: Path,
    reference_a: str,
    reference_b: str,
) -> tuple[Path, Path, str, str]:
    structure_a = pdb_ca_sequence(listed_a)
    structure_b = pdb_ca_sequence(listed_b)
    if structure_a == reference_a and structure_b == reference_b:
        return listed_a, listed_b, structure_a, structure_b
    if structure_a == reference_b and structure_b == reference_a:
        return listed_b, listed_a, structure_b, structure_a
    raise ValueError(
        f"Listed PDB sequences do not match either manual-row order for {case_id}: "
        f"{listed_a} ({len(structure_a)}), {listed_b} ({len(structure_b)}), "
        f"references ({len(reference_a)}, {len(reference_b)})"
    )

# ============================================================
# 构建 MALISAM / MALIDUP benchmark examples
# ============================================================

def mali_examples(dataset: str) -> list[dict]:
    rows = []
    root = MALI_ROOT / dataset
    pdb_pairs = mali_pdb_pairs(dataset)
    seen_cases = set()
    for manual_path in sorted(root.glob("*/*.manual.ali")):
        stem = Path(str(manual_path)[: -len(".manual.ali")])
        aligned_a, aligned_b = read_alignment(manual_path)
        seq_a = reference_sequence(aligned_a)
        seq_b = reference_sequence(aligned_b)
        case_id = stem.parent.name
        if case_id not in pdb_pairs:
            raise KeyError(f"{case_id} from {manual_path} is missing from {MALI_PDB_LISTS[dataset]}")
        pdb_a, pdb_b, struct_a, struct_b = order_listed_mali_pdbs(
            case_id, *pdb_pairs[case_id], seq_a, seq_b
        )
        map_a = sequence_index_map(struct_a, seq_a)
        map_b = sequence_index_map(struct_b, seq_b)
        seen_cases.add(case_id)
        rows.append({
            "benchmark": dataset,
            "id": case_id,
            "legacy_output_id": stem.name,
            "manual_path": str(manual_path),
            "truth": alignment_pairs(aligned_a.upper(), aligned_b.upper()),
            "uppercase_truth": alignment_pairs(aligned_a, aligned_b, uppercase_only=True),
            "seq_a": seq_a,
            "seq_b": seq_b,
            "pdb_a": pdb_a,
            "pdb_b": pdb_b,
            "struct_seq_a": struct_a,
            "struct_seq_b": struct_b,
            "map_a": map_a,
            "map_b": map_b,
        })
    unreferenced = sorted(set(pdb_pairs) - seen_cases)
    if unreferenced:
        raise ValueError(
            f"{MALI_PDB_LISTS[dataset]} contains cases without manual alignments: {unreferenced[:10]}"
        )
    return rows


# ============================================================
# 构建 HOMSTRAD examples
# ============================================================

def homstrad_examples(split_root: Path) -> list[dict]:
    rows = []
    for ali_path in sorted(HOMSTRAD_DATA.glob("*/*.ali")):
        family = ali_path.parent.name
        order, aligned = parse_homstrad_ali(ali_path)
        query, target = order[0], order[-1]
        references = {
            query: reference_sequence(aligned[query]),
            target: reference_sequence(aligned[target]),
        }
        selected_dir = split_root / family
        selected_dir.mkdir(parents=True, exist_ok=True)
        paths, chain_map = prepare_homstrad_structures(
            ali_path.parent / f"{family}-sup.pdb",
            references,
            selected_dir,
        )
        struct_a = pdb_ca_sequence(paths[query])
        struct_b = pdb_ca_sequence(paths[target])
        truth = alignment_pairs(aligned[query], aligned[target], uppercase_only=True)
        rows.append({
            "benchmark": "homstrad",
            "id": family,
            "query": query,
            "target": target,
            "query_chain": chain_map[query],
            "target_chain": chain_map[target],
            "truth": truth,
            "uppercase_truth": truth,
            "seq_a": references[query],
            "seq_b": references[target],
            "pdb_a": paths[query],
            "pdb_b": paths[target],
            "struct_seq_a": struct_a,
            "struct_seq_b": struct_b,
            "map_a": sequence_index_map(struct_a, references[query]),
            "map_b": sequence_index_map(struct_b, references[target]),
        })
    return rows


# ============================================================
# Main_function_model: finish a benchmark case：
# run DeepAlign -> parse alignment -> calculate F1
# ===========================================================
def evaluate_one(
    example: dict,
    binary: Path,
    alignment_root: Path,
    overwrite: bool,
    reuse_only: bool,
) -> dict:
    started = time.time() #time start
    benchmark = example["benchmark"]
    output_id = example.get("legacy_output_id", example["id"]) if reuse_only else example["id"]
    safe_id = output_id.replace("/", "_")
    result = {
        "benchmark": benchmark,
        "id": example["id"],
        "status": "ok",
        "pdb_a": str(example["pdb_a"]),
        "pdb_b": str(example["pdb_b"]),
        "len_a_reference": len(example["seq_a"]),
        "len_b_reference": len(example["seq_b"]),
        "len_a_structure": len(example["struct_seq_a"]),
        "len_b_structure": len(example["struct_seq_b"]),
        "mapped_a": len(example["map_a"]),
        "mapped_b": len(example["map_b"]),
        "truth_matches": len(example["truth"]),
        "uppercase_truth_matches": len(example["uppercase_truth"]),
    }
    try:
        run = run_deepalign(
            binary,
            example["pdb_a"],
            example["pdb_b"],
            safe_id,
            alignment_root / benchmark / "alignments",
            overwrite,
            reuse_only,
        )
        if run["status"] == "error":
            raise RuntimeError(run)
        predicted = predicted_pairs_from_fasta(run["fasta"], example["map_a"], example["map_b"])
        result.update(metrics(example["truth"], predicted))
        result.update({
            f"uppercase_{key}": value
            for key, value in metrics(example["uppercase_truth"], predicted).items()
        })
        result.update(parse_deepalign_score(run["score"]))
        result["fasta_path"] = str(run["fasta"])
        result["score_path"] = str(run["score"])
    except Exception as error:  
        result["status"] = "error"
        result["error"] = repr(error)
        for key, value in {
            "tp": 0,
            "fp": 0,
            "fn": 0,
            "precision": 0.0,
            "recall": 0.0,
            "sensitivity": 0.0,
            "f1": 0.0,
            "matches": 0,
        }.items():
            result[key] = value
            result[f"uppercase_{key}"] = value
    result["elapsed_seconds"] = time.time() - started
    return result


def describe(values: list[float]) -> dict:
    ordered = sorted(float(value) for value in values if math.isfinite(float(value)))
    if not ordered:
        return {"mean": float("nan"), "std": float("nan"), "median": float("nan")}
    return {
        "mean": float(mean(ordered)),
        "std": float(np.std(ordered, ddof=1)) if len(ordered) > 1 else 0.0,
        "median": float(median(ordered)),
        "min": ordered[0],
        "max": ordered[-1],
    }


def summarize_pair_metrics(ok_rows: list[dict], prefix: str = "") -> dict:
    tp = sum(row[f"{prefix}tp"] for row in ok_rows)
    fp = sum(row[f"{prefix}fp"] for row in ok_rows)
    fn = sum(row[f"{prefix}fn"] for row in ok_rows)
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    return {
        f"{prefix}macro_precision": describe([row[f"{prefix}precision"] for row in ok_rows]),
        f"{prefix}macro_recall": describe([row[f"{prefix}recall"] for row in ok_rows]),
        f"{prefix}macro_f1": describe([row[f"{prefix}f1"] for row in ok_rows]),
        f"{prefix}micro_precision": precision,
        f"{prefix}micro_recall": recall,
        f"{prefix}micro_f1": (
            2 * precision * recall / (precision + recall) if precision + recall else 0.0
        ),
        f"{prefix}bad_cases_f1_lt_0.5": int(
            sum(row[f"{prefix}f1"] < 0.5 for row in ok_rows)
        ),
    }


def summarize(rows: list[dict]) -> dict:

    ok_rows = [row for row in rows if row["status"] == "ok"]
    if len(ok_rows)!=len(rows):
        raise RuntimeError(f"run invalid") # 按道理不会存在failed status
    result = {
        "pairs": len(rows),
        "ok_pairs": len(ok_rows),
        "errors": len(rows) - len(ok_rows),
        "tm_score": describe([row.get("tm_score", float("nan")) for row in ok_rows]),
        "deep_score": describe([row.get("deep_score", float("nan")) for row in ok_rows]),
    }
    result.update(summarize_pair_metrics(ok_rows))
    result.update(summarize_pair_metrics(ok_rows, "uppercase_"))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmarks", nargs="+", default=["malisam", "malidup", "homstrad"])
    parser.add_argument("--binary", type=Path, default=DEFAULT_DEEPALIGN)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument(
        "--reuse-alignments-from",
        type=Path,
        help="Read existing per-case FASTA/score files without modifying them.",
    )
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--overwrite", action="store_true")
   
    args = parser.parse_args()

    args.binary = args.binary.resolve()
    binary_version = args.binary.name
    if args.output_root is None:
        args.output_root = REPO / "test_results" / f"full_alignment_{binary_version}"
    else:
        args.output_root = args.output_root.resolve()

    if args.reuse_alignments_from is not None:
        args.reuse_alignments_from = args.reuse_alignments_from.resolve()

    if not args.binary.is_file():
        raise FileNotFoundError(args.binary)

    alignment_root = args.reuse_alignments_from or args.output_root
    reuse_only = args.reuse_alignments_from is not None

    examples: list[dict] = []
    if "malisam" in args.benchmarks:
        examples.extend(mali_examples("malisam"))
    if "malidup" in args.benchmarks:
        examples.extend(mali_examples("malidup"))
    if "homstrad" in args.benchmarks:
        examples.extend(homstrad_examples(HOMSTRAD_CHAIN_Dir))
    if args.limit:
        examples = examples[: args.limit]
    if not examples:
        raise RuntimeError("No examples selected")

    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                evaluate_one,
                example,args.binary,alignment_root,args.overwrite,reuse_only,
        
            ): example
            for example in examples
        }
        for index, future in enumerate(as_completed(futures), 1):
            rows.append(future.result())
            if index % 50 == 0 or index == len(futures):
                ok = sum(row["status"] == "ok" for row in rows)
                print(f"[{index}/{len(futures)}] ok={ok} errors={len(rows)-ok}", flush=True)

    rows.sort(key=lambda row: (row["benchmark"], row["id"]))
    args.output_root.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_root / f"{binary_version}_alignment_benchmark.csv"
    json_path = args.output_root / f"{binary_version}_alignment_benchmark.json"
    fieldnames = sorted({key for row in rows for key in row})
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    groups = {"all": rows}
    for benchmark in sorted({row["benchmark"] for row in rows}):
        groups[benchmark] = [row for row in rows if row["benchmark"] == benchmark]
    report = {
        "deepalign": str(args.binary),
        "output_root": str(args.output_root),
        "alignment_root": str(alignment_root),
        "groups": {name: summarize(group_rows) for name, group_rows in groups.items()},
        "errors": [row for row in rows if row["status"] != "ok"][:100],
        "error_count": sum(row["status"] != "ok" for row in rows),
    }
    json_path.write_text(json.dumps(report, indent=2, allow_nan=True) + "\n")
    print(json.dumps(report, indent=2, allow_nan=True))


if __name__ == "__main__":
    main()
