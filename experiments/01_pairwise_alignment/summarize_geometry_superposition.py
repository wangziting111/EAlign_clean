#!/usr/bin/env python3
"""Summarize geometric quality from existing pairwise-alignment outputs.

No aligner is rerun.  The script reads saved EAlign, DeepAlign, TM-align,
DALI, GTalign, and Foldseek outputs, uses the evaluated uppercase-uppercase
manual columns as L_ref, and writes per-pair and aggregate audit tables.

EAlign, DeepAlign, TM-align, and GTalign use their native reported geometry.
DALI uses its reported transformation and equivalences.  Foldseek does not
report a rigid transformation, RMSD, or TM-score in the saved output, so these
quantities are computed after a least-squares Kabsch fit of its predicted
residue correspondences.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from statistics import mean

import numpy as np

import evaluate_pairwise_alignment as evaluation
from run_pairwise_alignment import BENCHMARKS, EXPECTED, load_cases, sha256


REPO = Path(__file__).resolve().parents[2]
PAIRWISE_ROOT = REPO / "test_results/pairwise_runs"
DEFAULT_OUTPUT = REPO / "test_results/01_pairwise_alignment/geometry"
METHODS = (
    "ealign",
    "deepalign",
    "tmalign",
    "dali",
    "gtalign_cuda",
    "gtalign_patched_cpu",
    "gtalign_patched_cpu_t20",
    "foldseek",
)
METHOD_LABELS = {
    "ealign": "EAlign",
    "deepalign": "DeepAlign",
    "tmalign": "TM-align",
    "dali": "DALI",
    "gtalign_cuda": "GTalign-CUDA",
    "gtalign_patched_cpu": "GTalign-CPU-1",
    "gtalign_patched_cpu_t20": "GTalign-CPU-20",
    "foldseek": "Foldseek",
}
GTALIGN_METHODS = {
    "gtalign_cuda",
    "gtalign_patched_cpu",
    "gtalign_patched_cpu_t20",
}

TMALIGN_SUMMARY = re.compile(
    r"Aligned length=\s*(\d+),\s*RMSD=\s*([0-9.eE+-]+)"
)
TMALIGN_LENGTH = re.compile(r"Length of Chain_([12]):\s*(\d+) residues")
TMALIGN_SCORE = re.compile(
    r"^TM-score=\s*([0-9.eE+-]+).*length of Chain_([12])", re.MULTILINE
)


def reference_length(row_a: str, row_b: str) -> int:
    """Count the uppercase-uppercase reference columns used by evaluation."""
    if len(row_a) != len(row_b):
        raise ValueError("Reference alignment rows have different lengths")
    length = sum(a.isupper() and b.isupper() for a, b in zip(row_a, row_b))
    if length == 0:
        raise ValueError("Reference alignment has no evaluated residue pairs")
    return length


def parse_deepalign_score(path: Path) -> tuple[int, float, float]:
    line = next(
        (
            raw.strip()
            for raw in path.read_text(errors="replace").splitlines()
            if raw.strip() and not raw.startswith("#")
        ),
        "",
    )
    blocks = [block.split() for block in line.split("->")]
    if len(blocks) < 3 or len(blocks[2]) < 3:
        raise ValueError(f"Malformed score file: {path}")
    return int(blocks[2][0]), float(blocks[2][1]), float(blocks[2][2])


def parse_tmalign(path: Path) -> tuple[int, float, float]:
    text = path.read_text(errors="replace")
    summary = TMALIGN_SUMMARY.search(text)
    lengths = {int(chain): int(length) for chain, length in TMALIGN_LENGTH.findall(text)}
    scores = {int(chain): float(score) for score, chain in TMALIGN_SCORE.findall(text)}
    if summary is None or set(lengths) != {1, 2} or set(scores) != {1, 2}:
        raise ValueError(f"Incomplete TM-align output: {path}")
    shorter = 1 if lengths[1] <= lengths[2] else 2
    return int(summary[1]), float(summary[2]), scores[shorter]


def parse_gtalign(path: Path) -> tuple[int, float, float]:
    payload = json.loads(path.read_text())
    search = payload.get("gtalign_search")
    hits = search.get("search_results", []) if isinstance(search, dict) else []
    if len(hits) != 1 or not isinstance(hits[0].get("hit_record"), dict):
        raise ValueError(f"Expected exactly one GTalign hit: {path}")
    hit = hits[0]["hit_record"]
    alignment = hit.get("alignment")
    if not isinstance(alignment, dict):
        raise ValueError(f"Missing GTalign alignment object: {path}")
    query_length = int(hit["query_length"])
    reference_length = int(hit["reference_length"])
    tm_score = (
        float(alignment["tmscore_query"])
        if query_length <= reference_length
        else float(alignment["tmscore_refn"])
    )
    return int(alignment["n_matched"]), float(alignment["rmsd"]), tm_score


def parse_dat_coordinates(path: Path) -> dict[int, tuple[float, float, float]]:
    coordinates = {}
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith("-acc"):
            continue
        fields = line.split()
        coordinates[int(fields[1])] = tuple(float(value) for value in fields[-3:])
    if not coordinates:
        raise ValueError(f"No -acc coordinates in {path}")
    return coordinates


def parse_dali_transform(path: Path) -> tuple[list[list[float]], list[float]]:
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("-matrix"):
            fields = line.replace('"', "").split()
            rows.append([float(value) for value in fields[-4:]])
    if len(rows) != 3:
        raise ValueError(f"Expected three DALI transformation rows: {path}")
    return [row[:3] for row in rows], [row[3] for row in rows]


def matrix_vector(
    matrix: list[list[float]], vector: tuple[float, float, float]
) -> tuple[float, float, float]:
    return tuple(sum(row[index] * vector[index] for index in range(3)) for row in matrix)


def tm_d0(length: int) -> float:
    if length <= 15:
        return 0.5
    return max(0.5, 1.24 * math.pow(length - 15, 1.0 / 3.0) - 1.8)


def kabsch_metrics(case, pairs: set[tuple[int, int]]) -> tuple[int, float, float]:
    if not pairs:
        raise ValueError(f"Foldseek alignment has no residue pairs: {case.case_id}")
    _, query_coordinates = evaluation.pdb_ca_coordinates(case.pdb_a)
    _, target_coordinates = evaluation.pdb_ca_coordinates(case.pdb_b)
    ordered = sorted(pairs)
    query = np.asarray([query_coordinates[i] for i, _ in ordered], dtype=float)
    target = np.asarray([target_coordinates[j] for _, j in ordered], dtype=float)
    query_centered = query - query.mean(axis=0)
    target_centered = target - target.mean(axis=0)
    left, _, right_t = np.linalg.svd(target_centered.T @ query_centered)
    rotation = left @ right_t
    if np.linalg.det(rotation) < 0:
        left[:, -1] *= -1
        rotation = left @ right_t
    squared_distances = np.sum(
        (target_centered @ rotation - query_centered) ** 2, axis=1
    )
    rmsd = math.sqrt(float(np.mean(squared_distances)))
    normalization_length = min(len(query_coordinates), len(target_coordinates))
    d0 = tm_d0(normalization_length)
    tm_score = float(
        np.sum(1.0 / (1.0 + squared_distances / (d0 * d0)))
    ) / normalization_length
    return len(ordered), rmsd, tm_score


def parse_foldseek(case, case_dir: Path) -> tuple[int, float, float]:
    pairs = evaluation.foldseek_prediction(
        case_dir / "foldseek.tsv", (case.pdb_a, case.pdb_b)
    )
    return kabsch_metrics(case, pairs)


def parse_dali(case_dir: Path, result: dict) -> tuple[int, float, float]:
    dat_files = list((case_dir / "DAT").glob("*.dat"))
    query_dat = next((path for path in dat_files if path.name.startswith("q001")), None)
    target_dat = next((path for path in dat_files if path.name.startswith("t001")), None)
    if query_dat is None or target_dat is None:
        raise FileNotFoundError(f"Missing DALI DAT files: {case_dir}")
    query = parse_dat_coordinates(query_dat)
    target = parse_dat_coordinates(target_dat)
    matrix, translation = parse_dali_transform(case_dir / result["report"])

    pairs = []
    with (case_dir / "alignment_pairs.tsv").open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pairs.append((int(row["query_dali_index"]), int(row["target_dali_index"])))
    if len(pairs) != int(result["lali"]):
        raise ValueError(f"DALI lali/pair-count mismatch: {case_dir}")

    squared_distances = []
    for query_index, target_index in pairs:
        rotated = matrix_vector(matrix, target[target_index])
        transformed = tuple(rotated[index] + translation[index] for index in range(3))
        squared_distances.append(
            sum(
                (transformed[index] - query[query_index][index]) ** 2
                for index in range(3)
            )
        )

    rmsd = math.sqrt(sum(squared_distances) / len(squared_distances))
    if abs(rmsd - float(result["rmsd"])) > 0.051:
        raise ValueError(
            f"Recomputed DALI RMSD disagrees with rounded report: {case_dir}: "
            f"{rmsd:.6f} vs {result['rmsd']}"
        )
    normalization_length = min(int(result["ca_len_a"]), int(result["ca_len_b"]))
    d0 = tm_d0(normalization_length)
    tm_score = sum(1.0 / (1.0 + squared / (d0 * d0)) for squared in squared_distances)
    tm_score /= normalization_length
    return len(pairs), rmsd, tm_score


def saved_result(
    method: str, benchmark: str, case_id: str, pairwise_root: Path
) -> tuple[Path, dict]:
    if method == "dali":
        case_dir = pairwise_root / "dali" / benchmark / "pairs" / case_id
    else:
        case_dir = pairwise_root / method / benchmark / case_id
    path = case_dir / "result.json"
    if not path.is_file():
        raise FileNotFoundError(path)
    return case_dir, json.loads(path.read_text())


def validate_saved_inputs(case, result: dict, case_dir: Path) -> None:
    saved_hashes = result.get("sha256", {})
    if not saved_hashes:
        return
    for name, path in (("pdb_a", case.pdb_a), ("pdb_b", case.pdb_b)):
        if saved_hashes.get(name) != sha256(path):
            raise ValueError(f"Input changed since run: {case_dir}: {path}")


def collect_rows(pairwise_root: Path, methods: list[str]) -> list[dict]:
    cases = load_cases(list(BENCHMARKS))
    rows = []
    for case in cases:
        ref_length = reference_length(case.row_a, case.row_b)
        for method in methods:
            case_dir, result = saved_result(
                method, case.benchmark, case.case_id, pairwise_root
            )
            expected_status = "hit" if method == "dali" else "ok"
            row = {
                "dataset": case.benchmark,
                "case_id": case.case_id,
                "method": method,
                "method_label": METHOD_LABELS[method],
                "status": result.get("status", ""),
                "reference_aligned_length": ref_length,
                "predicted_aligned_length": "",
                "aligned_length_ratio": "",
                "rmsd_angstrom": "",
                "tm_score": "",
                "geometry_source": "",
            }
            if result.get("status") == expected_status:
                if method in ("ealign", "deepalign"):
                    metrics = parse_deepalign_score(case_dir / f"{case.case_id}.score")
                    geometry_source = "native output"
                elif method == "tmalign":
                    metrics = parse_tmalign(case_dir / "stdout.txt")
                    geometry_source = "native output"
                elif method == "dali":
                    metrics = parse_dali(case_dir, result)
                    geometry_source = "native transform and equivalences"
                elif method in GTALIGN_METHODS:
                    validate_saved_inputs(case, result, case_dir)
                    metrics = parse_gtalign(case_dir / "gtalign.json")
                    geometry_source = "native output"
                elif method == "foldseek":
                    validate_saved_inputs(case, result, case_dir)
                    metrics = parse_foldseek(case, case_dir)
                    geometry_source = "Kabsch refit of predicted pairs"
                else:
                    raise ValueError(f"Unsupported method: {method}")
                aligned_length, rmsd, tm_score = metrics
                row.update(
                    {
                        "predicted_aligned_length": aligned_length,
                        "aligned_length_ratio": aligned_length / ref_length,
                        "rmsd_angstrom": rmsd,
                        "tm_score": tm_score,
                        "geometry_source": geometry_source,
                    }
                )
            rows.append(row)
    return rows


def summarize(rows: list[dict], methods: list[str]) -> list[dict]:
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["dataset"], row["method"])].append(row)
    summaries = []
    for dataset in BENCHMARKS:
        for method in methods:
            selected = grouped[(dataset, method)]
            measured = [row for row in selected if row["tm_score"] != ""]
            if len(selected) != EXPECTED[dataset]:
                raise ValueError(
                    f"Expected {EXPECTED[dataset]} {dataset}/{method} rows, got {len(selected)}"
                )
            summaries.append(
                {
                    "dataset": dataset,
                    "method": method,
                    "method_label": METHOD_LABELS[method],
                    "requested_pairs": len(selected),
                    "aligned_pairs": len(measured),
                    "coverage": len(measured) / len(selected),
                    "geometry_source": (
                        measured[0]["geometry_source"] if measured else ""
                    ),
                    "mean_aligned_length_ratio_hit_only": mean(
                        float(row["aligned_length_ratio"]) for row in measured
                    ),
                    "mean_rmsd_angstrom_hit_only": mean(
                        float(row["rmsd_angstrom"]) for row in measured
                    ),
                    "mean_tm_score_hit_only": mean(float(row["tm_score"]) for row in measured),
                    "mean_tm_score_all_pairs_zero_missing": sum(
                        float(row["tm_score"]) for row in measured
                    )
                    / len(selected),
                }
            )
    return summaries


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def print_summary(rows: list[dict]) -> None:
    print(
        "Dataset   Method      aligned/N  mean Lali/Lref  mean RMSD (A)  "
        "mean TM-score  TM all-pairs*"
    )
    for row in rows:
        print(
            f"{row['dataset']:<9} {row['method_label']:<11} "
            f"{row['aligned_pairs']:>4}/{row['requested_pairs']:<4} "
            f"{row['mean_aligned_length_ratio_hit_only']:>15.3f} "
            f"{row['mean_rmsd_angstrom_hit_only']:>14.3f} "
            f"{row['mean_tm_score_hit_only']:>14.3f} "
            f"{row['mean_tm_score_all_pairs_zero_missing']:>14.3f}"
        )
    print("* Missing/no-hit/error alignments contribute TM-score 0 only in this column.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--methods",
        nargs="+",
        choices=METHODS,
        default=list(METHODS),
        help="Methods to summarize; use the four supplementary methods for an incremental run.",
    )
    parser.add_argument(
        "--pairwise-root",
        type=Path,
        default=PAIRWISE_ROOT,
        help="Shared run root containing method/benchmark and dali/benchmark outputs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
    )
    args = parser.parse_args()
    pairwise_root = args.pairwise_root.resolve()
    rows = collect_rows(pairwise_root, args.methods)
    summaries = summarize(rows, args.methods)
    output = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else DEFAULT_OUTPUT.parent
        / (
            DEFAULT_OUTPUT.name
            if args.methods == list(METHODS)
            else "geometry_supplement"
        )
    )
    write_csv(output / "per_pair.csv", rows)
    write_csv(output / "summary.csv", summaries)
    print_summary(summaries)
    print(f"per-pair: {output / 'per_pair.csv'}")
    print(f"summary:  {output / 'summary.csv'}")


if __name__ == "__main__":
    main()
