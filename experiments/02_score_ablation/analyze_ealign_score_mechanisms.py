#!/usr/bin/env python3
"""Detect coherent correspondence shifts in rescued MALISAM failures.

A rescued failure satisfies

    uppercase_f1(ablation) < 0.5 <= uppercase_f1(Full).

For each rescued failure, the script finds reference pairs recovered by Full
EAlign but mapped to a different target residue by the ablation.  A coherent
shifted block is a maximal run of at least three pairs for which the query,
reference-target, and ablated-target indices are consecutive and the ablated
target maintains a constant offset from the reference target.

The analysis is exhaustive over rescued MALISAM failures.  Its primary
mechanism statistic is the fraction of restored wrong-target corrections that
belong to coherent shifted blocks.  It also selects the maximum-F1-rescue case
for each ablation as a deterministic illustration.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple


REPO = Path(__file__).resolve().parents[2]
DEFAULT_INPUTS = {
    "full": REPO
    / "test_results/full_alignment_ealign/ealign_alignment_benchmark.csv",
    "without_vector": REPO
    / "test_results/full_alignment_ealign_without_vector"
    / "ealign_without_vector_alignment_benchmark.csv",
    "without_3di": REPO
    / "test_results/full_alignment_ealign_without_3di"
    / "ealign_without_3di_alignment_benchmark.csv",
}
ALIGNMENT_DIRS = {
    "full": REPO / "test_results/full_alignment_ealign/malisam/alignments",
    "without_vector": REPO
    / "test_results/full_alignment_ealign_without_vector/malisam/alignments",
    "without_3di": REPO
    / "test_results/full_alignment_ealign_without_3di/malisam/alignments",
}


@dataclass(frozen=True)
class Block:
    a_start: int
    a_end: int
    reference_b_start: int
    reference_b_end: int
    ablated_b_start: int
    ablated_b_end: int
    offset: int
    length: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full-csv", type=Path, default=DEFAULT_INPUTS["full"])
    parser.add_argument(
        "--without-vector-csv",
        type=Path,
        default=DEFAULT_INPUTS["without_vector"],
    )
    parser.add_argument(
        "--without-3di-csv",
        type=Path,
        default=DEFAULT_INPUTS["without_3di"],
    )
    parser.add_argument("--bad-threshold", type=float, default=0.5)
    parser.add_argument("--min-block-length", type=int, default=3)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO / "test_results/ealign_score_mechanisms",
    )
    args = parser.parse_args()
    if args.min_block_length < 1:
        parser.error("--min-block-length must be at least 1")
    return args


def load_f1(path: Path) -> Dict[str, float]:
    if not path.is_file():
        raise SystemExit(f"benchmark CSV not found: {path}")
    result: Dict[str, float] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"benchmark", "id", "status", "uppercase_f1"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise SystemExit(f"{path} is missing columns: {sorted(missing)}")
        for row in reader:
            if row["benchmark"].strip().lower() != "malisam":
                continue
            if row["status"].strip().lower() != "ok":
                raise SystemExit(f"non-ok MALISAM row in {path}: {row['id']}")
            result[row["id"].strip()] = float(row["uppercase_f1"])
    return result


def read_alignment(path: Path) -> Tuple[str, str]:
    if not path.is_file():
        raise SystemExit(f"alignment not found: {path}")
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if any(line.startswith(">") for line in lines):
        sequences: List[str] = []
        current: List[str] = []
        for line in lines:
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current))
                    current = []
            else:
                current.append(line)
        if current:
            sequences.append("".join(current))
    else:
        sequences = [line for line in lines if not line.startswith(("#", "$"))]
    if len(sequences) != 2 or len(sequences[0]) != len(sequences[1]):
        raise SystemExit(f"expected a two-sequence alignment in {path}")
    return sequences[0], sequences[1]


def alignment_map(
    sequences: Tuple[str, str], *, uppercase_reference_only: bool = False
) -> Dict[int, int]:
    position_a = position_b = 0
    result: Dict[int, int] = {}
    for residue_a, residue_b in zip(*sequences):
        if residue_a != "-":
            position_a += 1
        if residue_b != "-":
            position_b += 1
        aligned = residue_a != "-" and residue_b != "-"
        if uppercase_reference_only:
            aligned = aligned and residue_a.isupper() and residue_b.isupper()
        if aligned:
            result[position_a] = position_b
    return result


def coherent_blocks(
    shifted_rows: Sequence[Tuple[int, int, int, int]], min_length: int
) -> List[Block]:
    runs: List[List[Tuple[int, int, int, int]]] = []
    for row in shifted_rows:
        if runs:
            previous = runs[-1][-1]
            indices_are_consecutive = all(
                row[index] == previous[index] + 1 for index in range(3)
            )
            if indices_are_consecutive and row[3] == previous[3]:
                runs[-1].append(row)
                continue
        runs.append([row])

    blocks = []
    for run in runs:
        if len(run) < min_length:
            continue
        first, last = run[0], run[-1]
        blocks.append(
            Block(
                a_start=first[0],
                a_end=last[0],
                reference_b_start=first[1],
                reference_b_end=last[1],
                ablated_b_start=first[2],
                ablated_b_end=last[2],
                offset=first[3],
                length=len(run),
            )
        )
    return blocks


def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    scores = {
        "full": load_f1(args.full_csv),
        "without_vector": load_f1(args.without_vector_csv),
        "without_3di": load_f1(args.without_3di_csv),
    }
    pair_ids = set(scores["full"])
    for label, rows in scores.items():
        if set(rows) != pair_ids:
            raise SystemExit(f"MALISAM pair IDs differ for {label}")

    case_rows: List[Dict[str, object]] = []
    block_rows: List[Dict[str, object]] = []
    for variant in ("without_vector", "without_3di"):
        rescued_ids = [
            pair_id
            for pair_id in sorted(pair_ids)
            if scores[variant][pair_id] < args.bad_threshold
            and scores["full"][pair_id] >= args.bad_threshold
        ]
        for pair_id in rescued_ids:
            reference = alignment_map(
                read_alignment(
                    REPO / "dataset/Malisam" / pair_id / f"{pair_id}.manual.ali"
                ),
                uppercase_reference_only=True,
            )
            full = alignment_map(
                read_alignment(ALIGNMENT_DIRS["full"] / f"{pair_id}.fasta")
            )
            ablation = alignment_map(
                read_alignment(ALIGNMENT_DIRS[variant] / f"{pair_id}.fasta")
            )

            corrected = []
            shifted = []
            dropped = []
            for position_a, reference_b in sorted(reference.items()):
                if full.get(position_a) != reference_b:
                    continue
                ablated_b = ablation.get(position_a)
                if ablated_b == reference_b:
                    continue
                corrected.append((position_a, reference_b))
                if ablated_b is None:
                    dropped.append((position_a, reference_b))
                else:
                    shifted.append(
                        (
                            position_a,
                            reference_b,
                            ablated_b,
                            ablated_b - reference_b,
                        )
                    )

            blocks = coherent_blocks(shifted, args.min_block_length)
            block_support = sum(block.length for block in blocks)
            case_rows.append(
                {
                    "variant": variant,
                    "pair_id": pair_id,
                    "full_uppercase_f1": scores["full"][pair_id],
                    "ablation_uppercase_f1": scores[variant][pair_id],
                    "rescue_gain": scores["full"][pair_id]
                    - scores[variant][pair_id],
                    "corrected_reference_pairs": len(corrected),
                    "wrong_target_pairs": len(shifted),
                    "dropped_pairs": len(dropped),
                    "coherent_blocks": len(blocks),
                    "wrong_target_pairs_in_blocks": block_support,
                    "coherent_shift_coverage": (
                        block_support / len(shifted) if shifted else 0.0
                    ),
                }
            )
            for block in blocks:
                block_rows.append(
                    {
                        "variant": variant,
                        "pair_id": pair_id,
                        **asdict(block),
                    }
                )

    summary_rows = []
    maximum_rows = []
    for variant in ("without_vector", "without_3di"):
        subset = [row for row in case_rows if row["variant"] == variant]
        wrong_target_pairs = sum(int(row["wrong_target_pairs"]) for row in subset)
        block_support = sum(
            int(row["wrong_target_pairs_in_blocks"]) for row in subset
        )
        summary_rows.append(
            {
                "variant": variant,
                "malisam_pairs": len(pair_ids),
                "rescued_failures": len(subset),
                "restored_reference_pairs": sum(
                    int(row["corrected_reference_pairs"]) for row in subset
                ),
                "wrong_target_corrections": wrong_target_pairs,
                "dropped_corrections": sum(
                    int(row["dropped_pairs"]) for row in subset
                ),
                "wrong_target_corrections_in_blocks": block_support,
                "coherent_shift_coverage": (
                    block_support / wrong_target_pairs if wrong_target_pairs else 0.0
                ),
            }
        )
        if subset:
            maximum_rows.append(
                min(
                    subset,
                    key=lambda row: (-float(row["rescue_gain"]), row["pair_id"]),
                )
            )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "rescued_case_metrics.csv", case_rows)
    write_csv(args.output_dir / "shift_blocks.csv", block_rows)
    write_csv(args.output_dir / "mechanism_summary.csv", summary_rows)
    write_csv(args.output_dir / "maximum_rescued_cases.csv", maximum_rows)
    (args.output_dir / "mechanism_summary.json").write_text(
        json.dumps(
            {
                "definitions": {
                    "reference_pairs": "uppercase-uppercase MALISAM columns",
                    "rescued_failure": (
                        f"uppercase_f1(ablation) < {args.bad_threshold} <= "
                        "uppercase_f1(Full)"
                    ),
                    "corrected_pair": (
                        "Full matches the reference and the ablation does not"
                    ),
                    "coherent_shifted_block": (
                        f"at least {args.min_block_length} consecutive corrected "
                        "wrong-target pairs with constant target offset"
                    ),
                    "coherent_shift_coverage": (
                        "restored wrong-target corrections in coherent shifted "
                        "blocks divided by all restored wrong-target corrections"
                    ),
                    "maximum_case": (
                        "largest uppercase_f1(Full)-uppercase_f1(ablation); "
                        "pair ID breaks exact ties"
                    ),
                },
                "summary": summary_rows,
                "maximum_rescued_cases": maximum_rows,
            },
            indent=2,
        )
        + "\n"
    )

    print("variant          rescued  restored  wrong-target  in-block  coverage  maximum rescue")
    for summary, maximum in zip(summary_rows, maximum_rows):
        print(
            f"{summary['variant']:<16} "
            f"{summary['rescued_failures']:>7} "
            f"{summary['restored_reference_pairs']:>9} "
            f"{summary['wrong_target_corrections']:>13} "
            f"{summary['wrong_target_corrections_in_blocks']:>9} "
            f"{summary['coherent_shift_coverage']:>9.3f}  "
            f"{maximum['pair_id']}"
        )
    print(f"\nWrote outputs to {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
