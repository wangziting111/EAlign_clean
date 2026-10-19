#!/usr/bin/env python3
"""Reconstruct and audit SABmark 1.65 from HF membership and SCOP/ASTRAL.

The script intentionally treats the Hugging Face files as pair-level metadata
only.  Coordinates and the independent ATOM-sequence reference come from the
historical SCOP/ASTRAL 1.65 release.
"""

import argparse
import csv
import hashlib
import json
import os
import shutil
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from difflib import SequenceMatcher
from pathlib import Path


EXPECTED = {
    "sup": {"groups": 425, "domains": 3280, "pairs": 19092},
    "twi": {"groups": 209, "domains": 1740, "pairs": 10667},
}

SOURCES = {
    "sup.jsonl": (
        "https://huggingface.co/datasets/DeepFoldProtein/SABmark-dataset/"
        "resolve/main/sup.jsonl?download=true",
        "93ebb3c713f727d6fd0b1d232c47d734b00f0c69dc66156bb14252748985d724",
        False,
    ),
    "twi.jsonl": (
        "https://huggingface.co/datasets/DeepFoldProtein/SABmark-dataset/"
        "resolve/main/twi.jsonl?download=true",
        "5be6df9996d1ef6abe26564945547fa0f5e1f37e66b2a49f7b66f4f6388e82b7",
        False,
    ),
    "dir.cla.scop.1.65.txt": (
        "https://scop.berkeley.edu/downloads/parse/dir.cla.scop.1.65.txt",
        "634b4e23009128669dda92ac8f6ef33ecc54fd6c02b3e71398c826f2b73feceb",
        True,
    ),
    "astral-scopdom-atom-gd-all-1.65.fa": (
        "https://scop.berkeley.edu/downloads/scopseq-1.65/"
        "astral-scopdom-atom-gd-all-1.65.fa",
        "3ab6d26c524520f94d11d1b7bc75dc0c36d419562e9945e40fd9f14bfca42676",
        True,
    ),
}

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "ASX": "B", "GLX": "Z", "XLE": "J", "SEC": "U", "PYL": "O",
    "UNK": "X", "MSE": "M", "FME": "M", "HYP": "P", "SEP": "S",
    "TPO": "T", "PTR": "Y", "CSO": "C", "CSD": "C", "CSE": "C",
    "CME": "C", "CYX": "C", "KCX": "K", "LLP": "K", "MLY": "K",
    "M3L": "K", "PCA": "E", "HIC": "H", "HIE": "H", "HID": "H",
    "HIP": "H", "ASH": "D", "GLH": "E", "AIB": "A", "DAL": "A",
    "DLE": "L", "DVA": "V", "SAR": "G",
}


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def urlopen(url, timeout=120, insecure=False):
    request = urllib.request.Request(url, headers={"User-Agent": "EAlign-SABmark/1.0"})
    context = ssl._create_unverified_context() if insecure else None
    return urllib.request.urlopen(request, timeout=timeout, context=context)


def download_atomic(url, destination, insecure=False, expected_sha=None, retries=4):
    destination = Path(destination)
    if destination.is_file():
        if expected_sha is None or sha256(destination) == expected_sha:
            return "existing"
        raise RuntimeError(f"checksum mismatch for existing file: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    last_error = None
    for attempt in range(1, retries + 1):
        temp = destination.with_name(destination.name + f".part.{os.getpid()}")
        try:
            with urlopen(url, insecure=insecure) as response, open(temp, "wb") as output:
                shutil.copyfileobj(response, output, length=1024 * 1024)
            if expected_sha and sha256(temp) != expected_sha:
                raise RuntimeError(f"checksum mismatch after downloading {url}")
            os.replace(temp, destination)
            return "downloaded"
        except Exception as exc:
            last_error = exc
            try:
                temp.unlink()
            except FileNotFoundError:
                pass
            if attempt < retries:
                time.sleep(min(2 ** attempt, 10))
    raise RuntimeError(f"failed to download {url}: {last_error}")


def fetch_sources(root):
    for name, (url, checksum, insecure) in SOURCES.items():
        status = download_atomic(url, root / name, insecure, checksum)
        print(f"{status:10s} {name}  sha256={checksum}")


def clean_sid(raw_sid):
    sid = raw_sid.split(":", 1)[-1]
    if len(sid) != 7 or not sid.startswith("d"):
        raise ValueError(f"unexpected SCOP SID: {raw_sid!r}")
    return sid


def load_hf_set(path, expected_name):
    groups = defaultdict(set)
    sequences = {}
    rows = []
    seen_pairs = set()
    with open(path, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            record = json.loads(line)
            if record.get("set_name") != expected_name:
                raise ValueError(f"{path}:{line_number}: wrong set_name")
            group = record["group_id"]
            sid1 = clean_sid(record["seq1_id"])
            sid2 = clean_sid(record["seq2_id"])
            if sid1 == sid2:
                raise ValueError(f"{path}:{line_number}: self pair {sid1}")
            canonical = (group,) + tuple(sorted((sid1, sid2)))
            if canonical in seen_pairs:
                raise ValueError(f"{path}:{line_number}: duplicate pair {canonical}")
            seen_pairs.add(canonical)
            for sid, sequence in ((sid1, record["seq1"]), (sid2, record["seq2"])):
                sequence = sequence.upper()
                if sid in sequences and sequences[sid] != sequence:
                    raise ValueError(f"{path}:{line_number}: conflicting sequence for {sid}")
                sequences[sid] = sequence
                groups[group].add(sid)
            rows.append({
                "set_name": expected_name,
                "group_id": group,
                "pair_id": record.get("pair_id", ""),
                "sid1": sid1,
                "sid2": sid2,
                "percent_identity": record.get("percent_identity", ""),
            })
    for group, sids in groups.items():
        wanted = len(sids) * (len(sids) - 1) // 2
        actual = sum(1 for row in rows if row["group_id"] == group)
        if actual != wanted:
            raise ValueError(f"{path}: group {group} has {actual} pairs; expected {wanted}")
    return groups, sequences, rows


def load_fasta(path):
    result = {}
    key = None
    parts = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if key is not None:
                    result[key] = "".join(parts).upper()
                key = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if key is not None:
        result[key] = "".join(parts).upper()
    return result


def load_scop_classification(path):
    result = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 6:
                result[fields[0]] = {
                    "pdb": fields[1], "region": fields[2], "sccs": fields[3],
                    "sunid": fields[4], "lineage": fields[5],
                }
    return result


def match_details(left, right):
    if left == right:
        return "exact", 1.0, 0
    if len(left) == len(right):
        differences = sum(a != b for a, b in zip(left, right))
        return "equal_length_mismatch", 1.0 - differences / max(1, len(left)), differences
    matcher = SequenceMatcher(None, left, right, autojunk=False)
    matches = sum(block.size for block in matcher.get_matching_blocks())
    if left in right:
        kind = "left_is_subsequence"
    elif right in left:
        kind = "right_is_subsequence"
    else:
        kind = "length_and_sequence_mismatch"
    return kind, matches / max(len(left), len(right), 1), max(len(left), len(right)) - matches


def write_tsv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def audit(root):
    manifests = root / "manifests"
    reports = root / "reports"
    manifests.mkdir(exist_ok=True)
    reports.mkdir(exist_ok=True)
    astral = load_fasta(root / "astral-scopdom-atom-gd-all-1.65.fa")
    scop = load_scop_classification(root / "dir.cla.scop.1.65.txt")
    summary = {"sets": {}, "cross_set": {}}
    all_sequences = {}
    set_sids = {}
    precheck_rows = []
    for set_name in ("sup", "twi"):
        groups, sequences, pairs = load_hf_set(root / f"{set_name}.jsonl", set_name)
        observed = {"groups": len(groups), "domains": len(sequences), "pairs": len(pairs)}
        if observed != EXPECTED[set_name]:
            raise RuntimeError(f"{set_name} audit failed: observed {observed}; expected {EXPECTED[set_name]}")
        set_sids[set_name] = set(sequences)
        domain_rows = []
        group_for_sid = {}
        for group, sids in groups.items():
            for sid in sids:
                if sid in group_for_sid:
                    raise RuntimeError(f"{set_name}: {sid} appears in multiple groups")
                group_for_sid[sid] = group
        for sid in sorted(sequences):
            if sid in all_sequences and all_sequences[sid] != sequences[sid]:
                raise RuntimeError(f"cross-set sequence conflict for {sid}")
            all_sequences[sid] = sequences[sid]
            ref = astral.get(sid, "")
            kind, identity, differences = match_details(sequences[sid], ref) if ref else ("missing", 0.0, len(sequences[sid]))
            row = {
                "set_name": set_name,
                "group_id": group_for_sid[sid],
                "sid": sid,
                "hf_length": len(sequences[sid]),
                "hf_sequence": sequences[sid],
                "scop_1_65_present": int(sid in scop),
                "astral_atom_present": int(sid in astral),
                "hf_vs_astral_atom": kind,
                "identity": f"{identity:.6f}",
                "differences": differences,
            }
            domain_rows.append(row)
            if kind != "exact" or sid not in scop:
                precheck_rows.append(row)
        write_tsv(
            manifests / f"{set_name}_domains.tsv",
            ["set_name", "group_id", "sid", "hf_length", "hf_sequence", "scop_1_65_present",
             "astral_atom_present", "hf_vs_astral_atom", "identity", "differences"],
            domain_rows,
        )
        write_tsv(
            manifests / f"{set_name}_pairs.tsv",
            ["set_name", "group_id", "pair_id", "sid1", "sid2", "percent_identity"],
            pairs,
        )
        summary["sets"][set_name] = {
            **observed,
            "expected": EXPECTED[set_name],
            "complete_pair_combinations": True,
            "missing_in_scop_1_65": sum(sid not in scop for sid in sequences),
            "missing_in_astral_atom_fasta": sum(sid not in astral for sid in sequences),
            "hf_astral_exact": sum(sequences[sid] == astral.get(sid) for sid in sequences),
            "hf_astral_mismatch": sum(sid in astral and sequences[sid] != astral[sid] for sid in sequences),
        }
    summary["cross_set"] = {
        "unique_domains": len(all_sequences),
        "overlap_domains": len(set_sids["sup"] & set_sids["twi"]),
        "total_reference_pairs": EXPECTED["sup"]["pairs"] + EXPECTED["twi"]["pairs"],
        "obsolete_or_unresolved_sids_in_scop_1_65": sum(sid not in scop for sid in all_sequences),
    }
    write_tsv(
        reports / "sequence_precheck.tsv",
        ["set_name", "group_id", "sid", "hf_length", "scop_1_65_present", "astral_atom_present",
         "hf_vs_astral_atom", "identity", "differences", "hf_sequence"],
        precheck_rows,
    )
    with open(reports / "membership_audit.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return all_sequences, astral, scop


def load_manifest_sequences(root):
    result = {}
    memberships = defaultdict(set)
    for set_name in ("sup", "twi"):
        path = root / "manifests" / f"{set_name}_domains.tsv"
        with open(path, encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                sid = row["sid"]
                if sid in result and result[sid] != row["hf_sequence"]:
                    raise RuntimeError(f"cross-set sequence conflict for {sid}")
                result[sid] = row["hf_sequence"]
                memberships[sid].add(set_name)
    return result, memberships


def structure_path(root, sid):
    return root / "structures" / sid[2:4] / f"{sid}.ent"


def validate_coordinate_file(path, sid):
    has_sid = False
    ca_count = 0
    with open(path, encoding="ascii", errors="replace") as handle:
        for line in handle:
            if sid in line[:100]:
                has_sid = True
            if line[:6].strip() in ("ATOM", "HETATM") and line[12:16].strip() == "CA":
                ca_count += 1
    if not has_sid:
        raise ValueError("downloaded file does not identify requested SID")
    if ca_count == 0:
        raise ValueError("downloaded file contains no CA coordinates")
    return ca_count


def coordinate_url(sid):
    quoted = urllib.parse.quote(sid, safe="._")
    return f"https://scop.berkeley.edu/downloads/pdbstyle/pdbstyle-1.65/{sid[2:4]}/{quoted}.ent"


def download_one_structure(root, sid, retries):
    destination = structure_path(root, sid)
    if destination.is_file():
        try:
            ca_count = validate_coordinate_file(destination, sid)
            return sid, "existing", ca_count, destination.stat().st_size, ""
        except Exception:
            destination.rename(destination.with_suffix(destination.suffix + ".invalid"))
    destination.parent.mkdir(parents=True, exist_ok=True)
    last_error = None
    for attempt in range(1, retries + 1):
        temp = destination.with_name(destination.name + f".part.{os.getpid()}.{os.getpid()}-{attempt}")
        try:
            with urlopen(coordinate_url(sid), insecure=True) as response, open(temp, "wb") as output:
                shutil.copyfileobj(response, output, length=1024 * 1024)
            ca_count = validate_coordinate_file(temp, sid)
            size = temp.stat().st_size
            os.replace(temp, destination)
            return sid, "downloaded", ca_count, size, ""
        except Exception as exc:
            last_error = exc
            try:
                temp.unlink()
            except FileNotFoundError:
                pass
            if attempt < retries:
                time.sleep(min(2 ** attempt, 10))
    return sid, "failed", 0, 0, str(last_error)


def download_structures(root, workers=8, retries=4):
    sequences, _ = load_manifest_sequences(root)
    sids = sorted(sequences)
    rows = []
    print(f"Downloading/verifying {len(sids)} unique SCOP 1.65 domain files with {workers} workers")
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(download_one_structure, root, sid, retries): sid for sid in sids}
        for index, future in enumerate(as_completed(futures), 1):
            sid, status, ca_count, size, error = future.result()
            rows.append({"sid": sid, "status": status, "ca_residues": ca_count, "bytes": size,
                         "path": str(structure_path(root, sid).relative_to(root)), "error": error})
            if index % 100 == 0 or index == len(sids):
                counts = Counter(row["status"] for row in rows)
                print(f"progress {index}/{len(sids)}  {dict(counts)}", flush=True)
    rows.sort(key=lambda row: row["sid"])
    write_tsv(root / "reports" / "structure_download.tsv",
              ["sid", "status", "ca_residues", "bytes", "path", "error"], rows)
    failed = [row for row in rows if row["status"] == "failed"]
    summary = {
        "requested_unique_domains": len(sids),
        "downloaded": sum(row["status"] == "downloaded" for row in rows),
        "existing_valid": sum(row["status"] == "existing" for row in rows),
        "failed": len(failed),
        "total_bytes": sum(row["bytes"] for row in rows),
    }
    with open(root / "reports" / "structure_download_summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    if failed:
        raise RuntimeError(f"{len(failed)} structure downloads failed; see reports/structure_download.tsv")


def coordinate_sequence(path):
    sequence = []
    seen_residues = set()
    unknown_names = Counter()
    saw_model = False
    with open(path, encoding="ascii", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip()
            if record == "MODEL":
                if saw_model:
                    break
                saw_model = True
                continue
            if record == "ENDMDL" and saw_model:
                break
            if record not in ("ATOM", "HETATM") or line[12:16].strip() != "CA":
                continue
            altloc = line[16:17]
            if altloc not in (" ", "A", "1"):
                continue
            residue_key = (line[21:22], line[22:26], line[26:27])
            if residue_key in seen_residues:
                continue
            seen_residues.add(residue_key)
            residue_name = line[17:20].strip().upper()
            amino_acid = AA3.get(residue_name)
            if amino_acid is None:
                amino_acid = "X"
                unknown_names[residue_name or "<blank>"] += 1
            sequence.append(amino_acid)
    return "".join(sequence), unknown_names


def verify_structures(root):
    hf_sequences, memberships = load_manifest_sequences(root)
    astral = load_fasta(root / "astral-scopdom-atom-gd-all-1.65.fa")
    scop = load_scop_classification(root / "dir.cla.scop.1.65.txt")
    rows = []
    for sid in sorted(hf_sequences):
        path = structure_path(root, sid)
        base = {
            "sid": sid,
            "sets": ",".join(sorted(memberships[sid])),
            "path": str(path.relative_to(root)),
            "scop_1_65_present": int(sid in scop),
            "file_present": int(path.is_file()),
            "hf_length": len(hf_sequences[sid]),
            "astral_atom_length": len(astral.get(sid, "")),
        }
        if not path.is_file():
            base.update({"coordinate_ca_length": 0, "hf_vs_coordinate": "missing",
                         "hf_coordinate_identity": "0.000000", "hf_coordinate_differences": len(hf_sequences[sid]),
                         "coordinate_vs_astral": "missing", "unknown_residue_names": ""})
            rows.append(base)
            continue
        coordinate, unknown = coordinate_sequence(path)
        hf_kind, hf_identity, hf_diff = match_details(hf_sequences[sid], coordinate)
        ar_kind, _, _ = match_details(coordinate, astral.get(sid, "")) if sid in astral else ("missing", 0.0, len(coordinate))
        base.update({
            "coordinate_ca_length": len(coordinate),
            "hf_vs_coordinate": hf_kind,
            "hf_coordinate_identity": f"{hf_identity:.6f}",
            "hf_coordinate_differences": hf_diff,
            "coordinate_vs_astral": ar_kind,
            "unknown_residue_names": ",".join(f"{key}:{value}" for key, value in sorted(unknown.items())),
        })
        rows.append(base)
    fields = ["sid", "sets", "path", "scop_1_65_present", "file_present", "hf_length",
              "astral_atom_length", "coordinate_ca_length", "hf_vs_coordinate", "hf_coordinate_identity",
              "hf_coordinate_differences", "coordinate_vs_astral", "unknown_residue_names"]
    write_tsv(root / "reports" / "structure_consistency.tsv", fields, rows)
    mismatch_rows = [
        row for row in rows
        if row["hf_vs_coordinate"] != "exact"
        or row["coordinate_vs_astral"] != "exact"
        or not row["file_present"]
        or not row["scop_1_65_present"]
        or bool(row["unknown_residue_names"])
    ]
    write_tsv(root / "reports" / "structure_exclusions_or_mismatches.tsv", fields, mismatch_rows)
    summary = {
        "unique_domains": len(rows),
        "files_present": sum(row["file_present"] for row in rows),
        "missing_files": sum(not row["file_present"] for row in rows),
        "exact_hf_vs_coordinate": sum(row["hf_vs_coordinate"] == "exact" for row in rows),
        "mismatch_hf_vs_coordinate": sum(row["file_present"] and row["hf_vs_coordinate"] != "exact" for row in rows),
        "coordinate_vs_astral_exact": sum(row["coordinate_vs_astral"] == "exact" for row in rows),
        "coordinate_vs_astral_mismatch": sum(row["file_present"] and row["coordinate_vs_astral"] != "exact" for row in rows),
        "unknown_coordinate_residue_names": sum(bool(row["unknown_residue_names"]) for row in rows),
        "obsolete_or_unresolved_sids_in_scop_1_65": sum(not row["scop_1_65_present"] for row in rows),
        "excluded_domains": sum(not row["file_present"] for row in rows),
        "domains_requiring_review": len(mismatch_rows),
    }
    with open(root / "reports" / "structure_consistency_summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    if summary["missing_files"]:
        raise RuntimeError("coordinate verification found missing files")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("fetch", "audit", "download", "verify", "all"))
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--retries", type=int, default=4)
    return parser.parse_args()


def main():
    args = parse_args()
    root = args.root.resolve()
    if args.command in ("fetch", "all"):
        fetch_sources(root)
    if args.command in ("audit", "all"):
        audit(root)
    if args.command in ("download", "all"):
        download_structures(root, max(1, args.workers), max(1, args.retries))
    if args.command in ("verify", "all"):
        verify_structures(root)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
