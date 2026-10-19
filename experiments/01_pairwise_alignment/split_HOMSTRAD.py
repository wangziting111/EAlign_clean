#!/usr/bin/env python3
"""Split the first and last .ali members from each HOMSTRAD *-sup.pdb."""

from __future__ import annotations

import argparse
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
HOMSTRAD = REPO / "data" / "HOMSTRAD"
DEFAULT_OUTPUT = REPO / "test_results/HOMSTRAD_split_structures"
AA3 = dict(zip(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL ASX GLX SEC PYL MSE".split(),
    "ARNDCQEGHILKMFPSTWYVBZUOM",
))
MISSING_LAST = {("cat", "1m85a"), ("hpr", "1ka5a")}


def parse_ali(path: Path) -> tuple[list[str], dict[str, str]]:
    order, rows = [], {}
    current = None
    skip_description = False
    for raw in path.read_text(errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("C;"):
            continue
        if line.startswith(">P1;"):
            current = line[4:].strip()
            if current in rows:
                raise ValueError(f"Duplicate member {current}: {path}")
            order.append(current)
            rows[current] = ""
            skip_description = True
        elif current is not None and skip_description:
            skip_description = False
        elif current is not None:
            rows[current] += line.replace("*", "")
    if len(order) < 2 or len({len(rows[name]) for name in order}) != 1:
        raise ValueError(f"Invalid HOMSTRAD alignment: {path}")
    return order, rows


def sequence(row: str) -> str:
    return "".join(c.upper() for c in row if c.isalpha())


def ca_sequence(path: Path, chain: str | None = None) -> str:
    residues = []
    previous = None
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or line[12:16].strip() != "CA":
            continue
        if chain is not None and line[21] != chain:
            continue
        aa = AA3.get(line[17:20].strip().upper())
        if aa is None:
            continue
        key = (line[21], line[22:26], line[26])
        if key != previous:
            residues.append(aa)
            previous = key
    return "".join(residues)


def chain_sequences(path: Path) -> dict[str, str]:
    chains = {line[21] for line in path.read_text(errors="replace").splitlines()
              if line.startswith(("ATOM  ", "HETATM")) and line[12:16].strip() == "CA"}
    return {chain: ca_sequence(path, chain) for chain in sorted(chains)}


def select_chain(path: Path, member: str, reference: str, chains: dict[str, str]) -> tuple[str, int]:
    family = path.parent.name
    matches = [c for c, seq in chains.items() if seq == reference]
    if not matches:
        matches = [c for c, seq in chains.items() if len(seq) > len(reference) and seq.startswith(reference)]
    if not matches and (family, member) in MISSING_LAST:
        matches = [c for c, seq in chains.items() if len(reference) == len(seq) + 1 and reference.startswith(seq)]
    if len(matches) != 1:
        raise ValueError(f"Expected one chain for {family}/{member}; found {matches} in {path}")
    chain = matches[0]
    return chain, min(len(reference), len(chains[chain]))


def write_trimmed_chain(sup: Path, output: Path, chain: str, length: int) -> None:
    lines = []
    previous = None
    index = -1
    for line in sup.read_text(errors="replace").splitlines(keepends=True):
        if not line.startswith(("ATOM  ", "HETATM")) or line[21] != chain:
            continue
        if line[17:20].strip().upper() not in AA3:
            continue
        key = (line[22:26], line[26])
        if key != previous:
            index += 1
            previous = key
        if index < length:
            lines.append(line)
    if index + 1 < length:
        raise ValueError(f"Only {index + 1}/{length} residues in {sup} chain {chain}")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("HEADER   \n" + "".join(lines) + "END\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--homstrad-root", type=Path, default=HOMSTRAD)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--limit", type=int, default=0, help="First N families for a smoke test")
    args = parser.parse_args()
    alignments = sorted(args.homstrad_root.glob("*/*.ali"))
    if not args.limit and len(alignments) != 1032:
        raise ValueError(f"Expected 1032 HOMSTRAD families, found {len(alignments)}")
    if args.limit:
        alignments = alignments[:args.limit]
    written = 0
    for ali in alignments:
        family = ali.parent.name
        order, rows = parse_ali(ali)
        members = (order[0], order[-1])
        sup = ali.parent / f"{family}-sup.pdb"
        chains = chain_sequences(sup)
        chosen = [select_chain(sup, member, sequence(rows[member]), chains) for member in members]
        if chosen[0][0] == chosen[1][0]:
            raise ValueError(f"Both members select the same chain: {family}")
        for member, (chain, length) in zip(members, chosen):
            target = args.output_root / family / f"{member}.pdb"
            expected = chains[chain][:length]
            write_trimmed_chain(sup, target, chain, length)
            if ca_sequence(target) != expected:
                raise ValueError(f"Written sequence mismatch: {target}")
            written += 1
    print(f"families={len(alignments)} structures_written={written}")


if __name__ == "__main__":
    main()
