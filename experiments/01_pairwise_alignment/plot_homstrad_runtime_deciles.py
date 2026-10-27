#!/usr/bin/env python3
"""Plot all 1,032 HOMSTRAD pair runtimes against L1*L2 with trend lines."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from run_pairwise_alignment import REPO, TIMING_SCOPE, ca_sequence, load_cases

CORE_METHODS = ("ealign", "deepalign", "tmalign")
SUPPLEMENTARY_METHODS = (
    "gtalign_cuda",
    "gtalign_patched_cpu",
    "gtalign_patched_cpu_t20",
    "foldseek",
)
ALL_METHODS = (*CORE_METHODS, *SUPPLEMENTARY_METHODS)
DEFAULT_ROOTS = {
    method: REPO / "test_results/pairwise_runs"
    for method in CORE_METHODS
}
DEFAULT_DALI_ROOT = REPO / "test_results/pairwise_runs/dali/homstrad"
DEFAULT_PAIRWISE_ROOT = REPO / "test_results/pairwise_runs"
DEFAULT_OUTPUT = REPO / "test_results/01_pairwise_alignment"
DEFAULT_DATA = DEFAULT_OUTPUT / "homstrad_runtime_scaling.csv"
DEFAULT_FIT = DEFAULT_OUTPUT / "homstrad_runtime_scaling_fit.csv"
DEFAULT_PNG = DEFAULT_OUTPUT / "homstrad_runtime_scaling.png"
DEFAULT_PDF = DEFAULT_OUTPUT / "homstrad_runtime_scaling.pdf"
DEFAULT_ALL_DATA = DEFAULT_OUTPUT / "homstrad_runtime_scaling_all_methods.csv"
DEFAULT_ALL_FIT = DEFAULT_OUTPUT / "homstrad_runtime_scaling_all_methods_fit.csv"
DEFAULT_ALL_PNG = DEFAULT_OUTPUT / "homstrad_runtime_scaling_all_methods.png"
DEFAULT_ALL_PDF = DEFAULT_OUTPUT / "homstrad_runtime_scaling_all_methods.pdf"

STYLES = {
    "ealign": ("EAlign", "#0072B2", "o"),
    "deepalign": ("DeepAlign", "#D55E00", "s"),
    "tmalign": ("TM-align", "#009E73", "^"),
    "dali": ("DALI", "#CC7987", "D"),
    "gtalign_cuda": ("GTalign-CUDA", "#7E57C2", "P"),
    "gtalign_patched_cpu": ("GTalign-CPU-1", "#E69F00", "v"),
    "gtalign_patched_cpu_t20": ("GTalign-CPU-20", "#56B4E9", "X"),
    "foldseek": ("Foldseek", "#4D4D4D", "h"),
}


def read_json(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(path)
    return json.loads(path.read_text())


def binary_runtime(root: Path, method: str, case_id: str) -> float:
    path = root / method / "homstrad" / case_id / "result.json"
    result = read_json(path)
    if result.get("method") != method or result.get("id") != case_id:
        raise ValueError(f"Metadata mismatch: {path}")
    if result.get("status") != "ok":
        raise ValueError(f"Unsuccessful {method} result: {path}")
    # The legacy core-method results predate the explicit timing_scope field,
    # but their `seconds` values use the same subprocess wall-clock boundary.
    if result.get("timing_scope") not in (None, TIMING_SCOPE):
        raise ValueError(f"Wrong timing scope: {path}")
    seconds = float(result["seconds"])
    if seconds <= 0:
        raise ValueError(f"Non-positive runtime: {path}")
    return seconds


def dali_runtime(root: Path, case_id: str) -> tuple[float, str]:
    path = root / "pairs" / case_id / "result.json"
    result = read_json(path)
    if result.get("family") != case_id:
        raise ValueError(f"Metadata mismatch: {path}")
    status = result.get("status")
    if status not in {"hit", "no_hit", "error", "timeout"}:
        raise ValueError(f"Unexpected Dali status {status!r}: {path}")
    seconds = float(result["seconds"])
    if seconds <= 0:
        raise ValueError(f"Non-positive runtime: {path}")
    return seconds, status


def collect_rows(
    split_root: Path,
    roots: dict[str, Path],
    dali_root: Path,
    methods: tuple[str, ...],
) -> list[dict]:
    cases = load_cases(["homstrad"], split_root)
    if len(cases) != 1032:
        raise ValueError(f"Expected 1032 HOMSTRAD pairs, found {len(cases)}")
    rows = []
    for case in cases:
        length_1 = len(ca_sequence(case.pdb_a))
        length_2 = len(ca_sequence(case.pdb_b))
        dali_seconds, dali_status = dali_runtime(dali_root, case.case_id)
        rows.append({
            "id": case.case_id,
            "L1": length_1,
            "L2": length_2,
            "L1_times_L2": length_1 * length_2,
            **{
                f"{method}_runtime_seconds": binary_runtime(roots[method], method, case.case_id)
                for method in methods
            },
            "dali_runtime_seconds": dali_seconds,
            "dali_status": dali_status,
        })
    return rows


def power_law_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Fit log10(t) = intercept + exponent*log10(L1*L2)."""
    log_x, log_y = np.log10(x), np.log10(y)
    exponent, intercept = np.polyfit(log_x, log_y, 1)
    predicted = intercept + exponent * log_x
    residual = np.sum((log_y - predicted) ** 2)
    total = np.sum((log_y - np.mean(log_y)) ** 2)
    r_squared = 1.0 - residual / total if total else 1.0
    return float(exponent), float(intercept), float(r_squared)


def write_csv(path: Path, rows: list[dict], fields: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({field: row[field] for field in fields} for row in rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ealign-root", type=Path, default=DEFAULT_ROOTS["ealign"])
    parser.add_argument("--deepalign-root", type=Path, default=DEFAULT_ROOTS["deepalign"])
    parser.add_argument("--tmalign-root", type=Path, default=DEFAULT_ROOTS["tmalign"])
    parser.add_argument("--dali-root", type=Path, default=DEFAULT_DALI_ROOT)
    parser.add_argument(
        "--all-methods",
        action="store_true",
        help="Also include GTalign-CUDA, GTalign-CPU-1/20, and Foldseek.",
    )
    parser.add_argument(
        "--pairwise-root",
        type=Path,
        default=DEFAULT_PAIRWISE_ROOT,
        help="Shared run root containing the supplementary method directories.",
    )
    parser.add_argument(
        "--split-root",
        type=Path,
        default=REPO / "test_results/HOMSTRAD_split_structures",
    )
    parser.add_argument("--data-output", type=Path)
    parser.add_argument("--fit-output", type=Path)
    parser.add_argument("--png", type=Path)
    parser.add_argument("--pdf", type=Path)
    args = parser.parse_args()

    methods = ALL_METHODS if args.all_methods else CORE_METHODS
    roots = {
        method: getattr(args, f"{method}_root").resolve()
        for method in CORE_METHODS
    }
    if args.all_methods:
        pairwise_root = args.pairwise_root.resolve()
        roots.update({method: pairwise_root for method in SUPPLEMENTARY_METHODS})

    data_output = args.data_output or (DEFAULT_ALL_DATA if args.all_methods else DEFAULT_DATA)
    fit_output = args.fit_output or (DEFAULT_ALL_FIT if args.all_methods else DEFAULT_FIT)
    png_output = args.png or (DEFAULT_ALL_PNG if args.all_methods else DEFAULT_PNG)
    pdf_output = args.pdf or (DEFAULT_ALL_PDF if args.all_methods else DEFAULT_PDF)

    rows = collect_rows(args.split_root, roots, args.dali_root.resolve(), methods)
    data_fields = [
        "id", "L1", "L2", "L1_times_L2",
        *(f"{method}_runtime_seconds" for method in methods),
        "dali_runtime_seconds", "dali_status",
    ]
    write_csv(data_output, rows, tuple(data_fields))

    figure_size = (8.5, 5.4) if args.all_methods else (7.0, 5.0)
    figure, axis = plt.subplots(figsize=figure_size, constrained_layout=True)
    x_all = np.asarray([row["L1_times_L2"] for row in rows], dtype=float)
    x_grid = np.geomspace(x_all.min(), x_all.max(), 300)
    fit_rows = []
    for method in (*methods, "dali"):
        label, color, marker = STYLES[method]
        valid = np.asarray([
            method != "dali" or row["dali_status"] in {"hit", "no_hit"}
            for row in rows
        ])
        x = x_all[valid]
        y = np.asarray([row[f"{method}_runtime_seconds"] for row in rows], dtype=float)[valid]
        exponent, intercept, r_squared = power_law_fit(x, y)
        axis.scatter(
            x, y, s=13, alpha=0.24, color=color, marker=marker,
            edgecolors="none", rasterized=True, label=f"{label} (n={len(x)})",
        )
        axis.plot(
            x_grid, (10 ** intercept) * x_grid ** exponent,
            color=color, linewidth=1.8,
        )
        fit_rows.append({
            "method": method,
            "points_in_trend": len(x),
            "power_law_exponent": exponent,
            "log10_intercept": intercept,
            "log10_r_squared": r_squared,
        })

    dali_failures = [row for row in rows if row["dali_status"] in {"error", "timeout"}]
    if dali_failures:
        axis.scatter(
            [row["L1_times_L2"] for row in dali_failures],
            [row["dali_runtime_seconds"] for row in dali_failures],
            s=28, marker="x", linewidths=1.0, color="#7A3E69",
            # label=f"Dali error/timeout (n={len(dali_failures)}; excluded from trend)",
            label=f"DALI parsing failure (n={len(dali_failures)})",
        )

    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel(r"Pair size product $L_1L_2$")
    axis.set_ylabel("Runtime (s)")
    title_suffix = " (all methods)" if args.all_methods else ""
    axis.set_title(f"Pairwise runtime scaling on HOMSTRAD{title_suffix}")
    axis.grid(True, which="major", linewidth=0.6, alpha=0.35)
    axis.grid(True, which="minor", linewidth=0.35, alpha=0.12)
    axis.legend(frameon=False, fontsize=8.2, ncol=3 if args.all_methods else 2)

    fit_fields = (
        "method", "points_in_trend", "power_law_exponent",
        "log10_intercept", "log10_r_squared",
    )
    write_csv(fit_output, fit_rows, fit_fields)
    for output in (png_output, pdf_output):
        output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(figure)

    print(f"pairs={len(rows)} dali_valid={len(rows) - len(dali_failures)} dali_failed={len(dali_failures)}")
    print(f"methods={','.join((*methods, 'dali'))}")
    print(f"data={data_output.resolve()}")
    print(f"fit={fit_output.resolve()}")
    print(f"png={png_output.resolve()}")
    print(f"pdf={pdf_output.resolve()}")


if __name__ == "__main__":
    main()
