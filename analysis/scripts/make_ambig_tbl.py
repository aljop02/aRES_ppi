#!/usr/bin/env python3
"""
make_ambig_tbl.py -- write a HADDOCK ambiguous interaction restraint (.tbl) file.

Exists because the HADDOCK 2.4 web server silently generated an EMPTY restraint
table for run `2_ares_T7` (job 780457): partner 1 was given 0 active residues and
304 passive, partner 2 was given 9 active, and the server produced no AIRs at all.
With `ranair=false`, `cmrest=false` and `surfrest=false`, nothing then drives the
partners together -- all 400 structures came back with the two chains ~270 A apart
and every intermolecular energy exactly zero. See docs/09 section 10.

The lesson: do not let the server infer the restraints. Write them, upload them,
and check the file is non-empty before spending queue time.

AIR semantics
-------------
One restraint per ACTIVE residue of the source molecule, targeting the union of
(active + passive) residues of the partner. HADDOCK evaluates each as a 1/r^6
summed effective distance over all atom pairs, so the standard `2.0 2.0 0.0`
bounds are satisfied as soon as the residue is in contact anywhere in the target
set. Passive residues are only ever targets, never sources.

Usage
-----
  # from a saved HADDOCK job parameter file (the reliable route)
  python3 make_ambig_tbl.py --from-params job_params_restrained.json -o ambig.tbl

  # or specify the selections directly
  python3 make_ambig_tbl.py --active B:103,104,105 --passive A:1-447 -o ambig.tbl

  # one-sided: only the named molecule's active residues become AIR sources
  python3 make_ambig_tbl.py --from-params params.json --source-segid B -o ambig.tbl
"""

from __future__ import annotations

import argparse
import json
import sys


def parse_selection(spec: str) -> tuple[str, list[int]]:
    """'B:103,104,111-114' -> ('B', [103, 104, 111, 112, 113, 114])"""
    if ":" not in spec:
        raise ValueError(f"selection must be SEGID:residues, got {spec!r}")
    segid, rest = spec.split(":", 1)
    out: list[int] = []
    for part in rest.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            lo, hi = part.split("-", 1)
            out.extend(range(int(lo), int(hi) + 1))
        else:
            out.append(int(part))
    return segid.strip(), sorted(set(out))


def load_params(path: str) -> dict[str, dict]:
    """Pull the active/passive lists and segids out of a HADDOCK job params file."""
    with open(path, encoding="utf-8") as fh:
        d = json.load(fh)
    mols = {}
    for key, p in d.get("partners", {}).items():
        segid = (p.get("segid") or p.get("chain") or key).strip()
        mols[segid] = {
            "active": sorted(set(p.get("activereslist") or [])),
            "passive": sorted(set(p.get("passivereslist") or [])),
        }
    return mols


def format_restraint(src_segid: str, resid: int,
                     tgt_segid: str, targets: list[int]) -> str:
    """One `assign` block: source residue -> the whole target selection."""
    lines = [f"assign ( resid {resid:>4} and segid {src_segid} )", "       ("]
    for i, t in enumerate(targets):
        joiner = "        " if i == 0 else "     or "
        lines.append(f"{joiner}( resid {t:>4} and segid {tgt_segid} )")
    lines.append("       )  2.0 2.0 0.0")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from-params", metavar="JSON",
                    help="HADDOCK job parameter file to read selections from")
    ap.add_argument("--active", action="append", default=[], metavar="SEGID:RESIDUES",
                    help="active residues, repeatable (e.g. B:103,104,111-114)")
    ap.add_argument("--passive", action="append", default=[], metavar="SEGID:RESIDUES",
                    help="passive residues, repeatable (e.g. A:1-447)")
    ap.add_argument("--source-segid", metavar="SEGID",
                    help="only this molecule's active residues become AIR sources "
                         "(one-sided docking); default is every molecule")
    ap.add_argument("-o", "--out", default="-", help="output .tbl (default stdout)")
    args = ap.parse_args()

    mols: dict[str, dict] = {}
    if args.from_params:
        mols = load_params(args.from_params)
    for spec in args.active:
        seg, res = parse_selection(spec)
        mols.setdefault(seg, {"active": [], "passive": []})["active"] = res
    for spec in args.passive:
        seg, res = parse_selection(spec)
        mols.setdefault(seg, {"active": [], "passive": []})["passive"] = res

    if len(mols) < 2:
        print(f"error: need >=2 molecules, found {len(mols)}", file=sys.stderr)
        return 2

    sources = [args.source_segid] if args.source_segid else list(mols)
    for s in sources:
        if s not in mols:
            print(f"error: segid {s!r} not among {list(mols)}", file=sys.stderr)
            return 2

    blocks: list[str] = []
    summary: list[str] = []
    for src in sources:
        active = mols[src]["active"]
        if not active:
            summary.append(f"  {src}: 0 active residues -> contributes NO restraints")
            continue
        for tgt, tm in mols.items():
            if tgt == src:
                continue
            targets = sorted(set(tm["active"]) | set(tm["passive"]))
            if not targets:
                print(f"error: molecule {tgt!r} has no active or passive residues, "
                      f"so restraints from {src!r} would have no target", file=sys.stderr)
                return 2
            for resid in active:
                blocks.append(format_restraint(src, resid, tgt, targets))
            summary.append(f"  {src} ({len(active)} active) -> {tgt} "
                           f"({len(targets)} target residues): {len(active)} restraints")

    if not blocks:
        print("error: produced 0 restraints -- this is exactly the failure that "
              "broke job 780457. Check that the source molecule has active residues.",
              file=sys.stderr)
        return 1

    header = [
        "! Ambiguous interaction restraints (AIRs)",
        "! Generated by analysis/scripts/make_ambig_tbl.py",
        f"! {len(blocks)} restraints",
    ]
    for line in summary:
        header.append(f"!{line}")
    text = "\n".join(header) + "\n!\n" + "\n".join(blocks) + "\n"

    if args.out == "-":
        sys.stdout.write(text)
    else:
        with open(args.out, "w", encoding="utf-8") as fh:
            fh.write(text)
        print(f"wrote {args.out}: {len(blocks)} restraints", file=sys.stderr)
    for line in summary:
        print(line, file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
