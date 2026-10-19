#!/usr/bin/env python3
"""Download SCOP 1.65 coordinate files needed only by SABmark TwiFP.

The script intentionally performs no sequence, chain/range, residue, or CA
validation. Existing coordinate files are left untouched.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import ssl
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


HF_URL = (
    "https://huggingface.co/datasets/DeepFoldProtein/SABmark-dataset/"
    "resolve/main/twi_fp.jsonl?download=true"
)
SCOP_URL = "https://scop.berkeley.edu/downloads/pdbstyle/pdbstyle-1.65"


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--twifp-jsonl", type=Path, default=here / "twi_fp.jsonl")
    parser.add_argument("--structures", type=Path, default=here.parent / "structures")
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--retries", type=int, default=4)
    return parser.parse_args()


def opener(url: str, timeout: int = 120):
    request = urllib.request.Request(url, headers={"User-Agent": "EAlign-SABmark/1.0"})
    return urllib.request.urlopen(
        request, timeout=timeout, context=ssl._create_unverified_context()
    )


def download_jsonl(path: Path, retries: int) -> None:
    if path.is_file():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        temporary = path.with_name(f"{path.name}.part.{os.getpid()}")
        try:
            with opener(HF_URL) as response, temporary.open("wb") as output:
                shutil.copyfileobj(response, output, length=1024 * 1024)
            os.replace(temporary, path)
            return
        except Exception as exc:
            last_error = exc
            temporary.unlink(missing_ok=True)
            if attempt < retries:
                time.sleep(min(2**attempt, 10))
    raise RuntimeError(f"cannot download {HF_URL}: {last_error}")


def clean_sid(value: str) -> str:
    sid = value.rsplit(":", 1)[-1]
    if len(sid) != 7 or not sid.startswith("d"):
        raise ValueError(f"unexpected SCOP SID: {value!r}")
    return sid


def load_sids(path: Path) -> set[str]:
    sids: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            record = json.loads(line)
            try:
                sids.add(clean_sid(record["seq1_id"]))
                sids.add(clean_sid(record["seq2_id"]))
            except (KeyError, ValueError) as exc:
                raise RuntimeError(f"{path}:{line_number}: {exc}") from exc
    return sids


def destination_for(structures: Path, sid: str) -> Path:
    return structures / sid[2:4] / f"{sid}.ent"


def url_for(sid: str) -> str:
    quoted = urllib.parse.quote(sid, safe="._")
    return f"{SCOP_URL}/{sid[2:4]}/{quoted}.ent"


def download_one(structures: Path, sid: str, retries: int) -> dict[str, object]:
    destination = destination_for(structures, sid)
    if destination.is_file():
        return {
            "sid": sid,
            "status": "existing",
            "bytes": destination.stat().st_size,
            "path": str(destination),
            "error": "",
        }

    destination.parent.mkdir(parents=True, exist_ok=True)
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        temporary = destination.with_name(
            f"{destination.name}.part.{os.getpid()}.{attempt}"
        )
        try:
            with opener(url_for(sid)) as response, temporary.open("wb") as output:
                shutil.copyfileobj(response, output, length=1024 * 1024)
            size = temporary.stat().st_size
            os.replace(temporary, destination)
            return {
                "sid": sid,
                "status": "downloaded",
                "bytes": size,
                "path": str(destination),
                "error": "",
            }
        except Exception as exc:
            last_error = exc
            temporary.unlink(missing_ok=True)
            if attempt < retries:
                time.sleep(min(2**attempt, 10))

    return {
        "sid": sid,
        "status": "failed",
        "bytes": 0,
        "path": str(destination),
        "error": str(last_error),
    }


def write_lines(path: Path, values: list[str]) -> None:
    path.write_text("".join(f"{value}\n" for value in values), encoding="utf-8")


def main() -> int:
    args = parse_args()
    twifp_jsonl = args.twifp_jsonl.resolve()
    structures = args.structures.resolve()
    workers = max(1, args.workers)
    retries = max(1, args.retries)

    download_jsonl(twifp_jsonl, retries)
    sids = sorted(load_sids(twifp_jsonl))
    missing = [sid for sid in sids if not destination_for(structures, sid).is_file()]

    report_dir = twifp_jsonl.parent / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    write_lines(report_dir / "twifp_domains.txt", sids)
    write_lines(report_dir / "twifp_missing_before_download.txt", missing)

    print(f"TwiFP unique domains: {len(sids)}")
    print(f"Already present:      {len(sids) - len(missing)}")
    print(f"To download:          {len(missing)}")

    rows: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(download_one, structures, sid, retries): sid
            for sid in missing
        }
        for index, future in enumerate(as_completed(futures), 1):
            rows.append(future.result())
            if index % 50 == 0 or index == len(missing):
                counts = Counter(str(row["status"]) for row in rows)
                print(f"progress {index}/{len(missing)} {dict(counts)}", flush=True)

    rows.sort(key=lambda row: str(row["sid"]))
    report_path = report_dir / "twifp_structure_download.tsv"
    with report_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("sid", "status", "bytes", "path", "error"),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    counts = Counter(str(row["status"]) for row in rows)
    failed = counts.get("failed", 0)
    summary = {
        "twifp_unique_domains": len(sids),
        "already_present_before_run": len(sids) - len(missing),
        "requested_downloads": len(missing),
        "downloaded": counts.get("downloaded", 0),
        "failed": failed,
        "validation_performed": False,
    }
    summary_path = report_dir / "twifp_structure_download_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    if failed:
        print(f"ERROR: {failed} downloads failed; see {report_path}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
