#!/usr/bin/env python3
"""Plot reference, ablated, and Full-EAlign patterns for maximum rescued cases.

The released MALISAM alignment defines the horizontal coordinate system.  Each
case is drawn as one correspondence diagram with three overlaid line sets:
released reference, ablation, and Full EAlign.  Residue dots are replaced by
the original reference-sequence letters.  Residues participating in the
released uppercase--uppercase reference correspondence are bold, while the
remaining residues are regular gray.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Mapping, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


REPO = Path(__file__).resolve().parents[2]
DEFAULT_CASES = REPO / "test_results/ealign_score_mechanisms/maximum_rescued_cases.csv"
FULL_DIR = REPO / "test_results/full_alignment_ealign/malisam/alignments"
VARIANT_DIRS = {
    "without_vector": REPO
    / "test_results/full_alignment_ealign_without_vector/malisam/alignments",
    "without_3di": REPO
    / "test_results/full_alignment_ealign_without_3di/malisam/alignments",
}
VARIANT_LABELS = {
    "without_vector": r"w/o orientation $V$",
    "without_3di": "w/o 3Di",
}

COLOR_ABLATION = "#D55E00"
COLOR_FULL = "#0072B2"
COLOR_REFERENCE = "#777777"
COLOR_MUTED = "#8A8A8A"
COLOR_TEXT = "#202020"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases-csv", type=Path, default=DEFAULT_CASES)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO / "test_results/ealign_score_mechanisms/figures",
    )
    parser.add_argument("--dpi", type=int, default=240)
    return parser.parse_args()


def read_alignment(path: Path) -> Tuple[List[str], Tuple[str, str]]:
    if not path.is_file():
        raise SystemExit(f"alignment not found: {path}")
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    headers: List[str] = []
    if any(line.startswith(">") for line in lines):
        sequences: List[str] = []
        current: List[str] = []
        for line in lines:
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current))
                    current = []
                headers.append(line[1:].strip())
            else:
                current.append(line)
        if current:
            sequences.append("".join(current))
    else:
        sequences = [line for line in lines if not line.startswith(("#", "$"))]
    if len(sequences) != 2 or len(sequences[0]) != len(sequences[1]):
        raise SystemExit(f"expected a two-sequence alignment in {path}")
    return headers, (sequences[0], sequences[1])


def alignment_map(sequences: Tuple[str, str]) -> Dict[int, int]:
    position_a = position_b = 0
    mapping: Dict[int, int] = {}
    for residue_a, residue_b in zip(*sequences):
        if residue_a != "-":
            position_a += 1
        if residue_b != "-":
            position_b += 1
        if residue_a != "-" and residue_b != "-":
            mapping[position_a] = position_b
    return mapping


def reference_layout(
    sequences: Tuple[str, str]
) -> Tuple[
    Dict[int, float],
    Dict[int, float],
    Dict[int, str],
    Dict[int, str],
    List[Tuple[int, int]],
]:
    position_a = position_b = 0
    x_a: Dict[int, float] = {}
    x_b: Dict[int, float] = {}
    letters_a: Dict[int, str] = {}
    letters_b: Dict[int, str] = {}
    uppercase_pairs: List[Tuple[int, int]] = []
    for column, (residue_a, residue_b) in enumerate(zip(*sequences), start=1):
        if residue_a != "-":
            position_a += 1
            x_a[position_a] = float(column)
            letters_a[position_a] = residue_a
        if residue_b != "-":
            position_b += 1
            x_b[position_b] = float(column)
            letters_b[position_b] = residue_b
        if residue_a.isupper() and residue_b.isupper():
            uppercase_pairs.append((position_a, position_b))
    return x_a, x_b, letters_a, letters_b, uppercase_pairs


def read_rows(path: Path) -> List[Dict[str, str]]:
    if not path.is_file():
        raise SystemExit(f"CSV not found: {path}")
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def draw_overlay(
    ax: plt.Axes,
    reference_map: Mapping[int, int],
    ablation_map: Mapping[int, int],
    full_map: Mapping[int, int],
    x_a: Mapping[int, float],
    x_b: Mapping[int, float],
    letters_a: Mapping[int, str],
    letters_b: Mapping[int, str],
    labels: Tuple[str, str],
) -> None:
    # Wide-to-narrow rendering leaves all three colors visible when mappings
    # coincide: gray reference halo, orange ablation, then blue Full EAlign.
    styles = (
        (reference_map, COLOR_REFERENCE, 3.0, 0.30, 1),
        (ablation_map, COLOR_ABLATION, 1.75, 0.76, 2),
        (full_map, COLOR_FULL, 0.75, 0.94, 3),
    )
    for mapping, color, linewidth, alpha, zorder in styles:
        for position_a, position_b in sorted(mapping.items()):
            if position_a not in x_a or position_b not in x_b:
                continue
            ax.plot(
                [x_a[position_a], x_b[position_b]],
                [0.90, 0.10],
                color=color,
                alpha=alpha,
                linewidth=linewidth,
                zorder=zorder,
            )

    aligned_a = set(reference_map)
    aligned_b = set(reference_map.values())
    for position, x in x_a.items():
        aligned = position in aligned_a
        ax.text(
            x,
            1.0,
            letters_a[position],
            ha="center",
            va="center",
            family="DejaVu Sans Mono",
            fontsize=7.0,
            fontweight="bold" if aligned else "normal",
            color=COLOR_TEXT if aligned else COLOR_MUTED,
            zorder=4,
        )
    for position, x in x_b.items():
        aligned = position in aligned_b
        ax.text(
            x,
            0.0,
            letters_b[position],
            ha="center",
            va="center",
            family="DejaVu Sans Mono",
            fontsize=7.0,
            fontweight="bold" if aligned else "normal",
            color=COLOR_TEXT if aligned else COLOR_MUTED,
            zorder=4,
        )
    ax.text(
        -0.012,
        1.0,
        labels[0],
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="center",
        fontsize=8,
        color=COLOR_TEXT,
    )
    ax.text(
        -0.012,
        0.0,
        labels[1],
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="center",
        fontsize=8,
        color=COLOR_TEXT,
    )
    ax.set_ylim(-0.25, 1.25)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_color(COLOR_REFERENCE)
    ax.tick_params(axis="x", colors=COLOR_MUTED, labelsize=7, length=2)


def main() -> None:
    args = parse_args()
    cases = read_rows(args.cases_csv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for case in cases:
        variant = case["variant"]
        pair_id = case["pair_id"]
        reference_path = REPO / "dataset/Malisam" / pair_id / f"{pair_id}.manual.ali"
        ablation_path = VARIANT_DIRS[variant] / f"{pair_id}.fasta"
        full_path = FULL_DIR / f"{pair_id}.fasta"

        _, reference = read_alignment(reference_path)
        ablation_headers, ablation = read_alignment(ablation_path)
        full_headers, full = read_alignment(full_path)
        headers = full_headers or ablation_headers or ["A", "B"]
        labels = (headers[0], headers[1])
        x_a, x_b, letters_a, letters_b, reference_pairs = reference_layout(reference)
        reference_map = dict(reference_pairs)
        ablation_map = alignment_map(ablation)
        full_map = alignment_map(full)

        columns = len(reference[0])
        width = max(11.5, columns * 0.125)
        fig, ax = plt.subplots(figsize=(width, 3.7))
        fig.subplots_adjust(left=0.105, right=0.99, top=0.80, bottom=0.18)
        fig.suptitle(
            (
                f"{pair_id}: reference, {VARIANT_LABELS[variant]}, and "
                "Full EAlign correspondences"
            ),
            x=0.105,
            ha="left",
            fontsize=11,
            fontweight="normal",
            color=COLOR_TEXT,
        )
        ax.set_title(
            (
                f"Ablation $F_1$={float(case['ablation_uppercase_f1']):.3f}; "
                f"Full EAlign $F_1$={float(case['full_uppercase_f1']):.3f}"
            ),
            loc="left",
            fontsize=9,
            fontweight="normal",
            color=COLOR_TEXT,
        )
        ax.set_xlim(0.25, columns + 0.75)
        draw_overlay(
            ax,
            reference_map,
            ablation_map,
            full_map,
            x_a,
            x_b,
            letters_a,
            letters_b,
            labels,
        )
        ax.set_xlabel("Released-reference alignment column", fontsize=8, color=COLOR_TEXT)
        tick_step = 10 if columns <= 120 else 20
        ax.set_xticks(range(tick_step, columns + 1, tick_step))
        legend = [
            Line2D([0], [0], color=COLOR_REFERENCE, alpha=0.55, lw=3.0, label="Reference"),
            Line2D([0], [0], color=COLOR_ABLATION, lw=1.75, label="Ablation"),
            Line2D([0], [0], color=COLOR_FULL, lw=0.9, label="Full EAlign"),
        ]
        fig.legend(
            handles=legend,
            loc="upper right",
            bbox_to_anchor=(0.99, 0.985),
            ncol=3,
            frameon=False,
            fontsize=8,
        )

        stem = args.output_dir / f"{pair_id}_alignment_patterns"
        fig.savefig(stem.with_suffix(".png"), dpi=args.dpi, bbox_inches="tight")
        fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
        fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
        plt.close(fig)
        print(stem.with_suffix(".png"))


if __name__ == "__main__":
    main()
