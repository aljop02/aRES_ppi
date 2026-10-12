#!/usr/bin/env python3
"""
zn_site_check.py -- find candidate Zn(Cys)4 / Zn(Cys,His) sites in a structure.

AlphaFold models carry no metals, but they place the coordinating side chains in
metal-binding geometry anyway. A tight cluster of four cysteine sulfurs is therefore
good evidence for a structural zinc site even with nothing modelled in the pocket.

The test is quantitative, not a guess. In an ideal tetrahedral Zn(Cys)4 site with a
Zn-S bond of 2.3 A, the S...S edge is

    2.3 * sqrt(8/3) = 3.76 A

so a 4-cysteine cluster with all pairwise SG-SG distances near 3.7-3.9 A is a zinc
site. Disulfides are a different thing entirely and much shorter: an S-S bond is
2.05 A. The script reports both so they cannot be confused.

Also prints mean pLDDT (B-factor) per region when --regions is given, which is how
you decide whether a loop or tail is worth docking at all.

Dependencies: numpy.

Usage
-----
  python3 zn_site_check.py model.pdb
  python3 zn_site_check.py model.pdb --regions 1-48:Zn-ribbon 49-120:NTD1 448-459:tail
"""

from __future__ import annotations

import argparse
import itertools
import sys
from collections import defaultdict

import numpy as np

IDEAL_ZN_SS = 2.3 * np.sqrt(8.0 / 3.0)  # 3.76 A
DISULFIDE_SS = 2.05


def load(path: str):
    """-> (coordinating atoms, per-residue mean B-factor, residue names)"""
    coord, bfac, names = {}, defaultdict(list), {}
    with open(path, errors="replace") as fh:
        for ln in fh:
            if not ln.startswith(("ATOM  ", "HETATM")):
                if ln.startswith("ENDMDL"):
                    break
                continue
            try:
                resseq = int(ln[22:26])
                xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            except ValueError:
                continue
            resname, atom, chain = ln[17:20].strip().upper(), ln[12:16].strip(), ln[21]
            try:
                bfac[(chain, resseq)].append(float(ln[60:66]))
            except ValueError:
                pass
            names[(chain, resseq)] = resname
            # Side-chain atoms that coordinate zinc.
            if (resname, atom) in (("CYS", "SG"), ("HIS", "NE2"), ("HIS", "ND1"),
                                   ("ASP", "OD1"), ("GLU", "OE1")):
                coord.setdefault((chain, resseq, resname), {})[atom] = xyz
    return coord, {k: float(np.mean(v)) for k, v in bfac.items()}, names


def main() -> None:
    ap = argparse.ArgumentParser(description="Detect candidate structural Zn sites.")
    ap.add_argument("structure", help="PDB file")
    ap.add_argument("--max-ss", type=float, default=4.5,
                    help="max pairwise SG-SG distance within a cluster (default 4.5 A)")
    ap.add_argument("--regions", nargs="*", default=[], metavar="LO-HI:NAME",
                    help="regions to report mean pLDDT for")
    args = ap.parse_args()

    coord, bfac, names = load(args.structure)
    cys = {k: v["SG"] for k, v in coord.items() if k[2] == "CYS" and "SG" in v}

    print(f"\n{args.structure}")
    print(f"  cysteines with SG: {len(cys)}")
    if not cys:
        print("  no cysteines -- nothing to test")
    else:
        keys = sorted(cys, key=lambda k: (k[0], k[1]))
        print(f"  positions: {', '.join(f'{k[0]}:{k[1]}' for k in keys)}")

        pairs = []
        for a, b in itertools.combinations(keys, 2):
            d = float(np.linalg.norm(cys[a] - cys[b]))
            if d < 9.0:
                pairs.append((a, b, d))
        if pairs:
            print(f"\n  SG-SG distances < 9 A:")
            for a, b, d in sorted(pairs, key=lambda x: x[2]):
                if d < 2.5:
                    tag = "  <- DISULFIDE"
                elif abs(d - IDEAL_ZN_SS) < 0.6:
                    tag = "  <- Zn-like"
                else:
                    tag = ""
                print(f"    {a[0]}:C{a[1]:<5} {b[0]}:C{b[1]:<5} {d:5.2f} A{tag}")

        print(f"\n  candidate Zn(Cys)4 sites (all 6 pairs < {args.max_ss} A):")
        found = 0
        for combo in itertools.combinations(keys, 4):
            ds = [float(np.linalg.norm(cys[a] - cys[b]))
                  for a, b in itertools.combinations(combo, 2)]
            if max(ds) < args.max_ss:
                found += 1
                label = ", ".join(f"{k[0]}:C{k[1]}" for k in combo)
                dev = abs(np.mean(ds) - IDEAL_ZN_SS)
                print(f"    {label}")
                print(f"       SG-SG range {min(ds):.2f}-{max(ds):.2f} A, "
                      f"mean {np.mean(ds):.2f} A "
                      f"(ideal tetrahedral Zn(Cys)4: {IDEAL_ZN_SS:.2f} A, "
                      f"deviation {dev:.2f} A)")
                spacing = [combo[i + 1][1] - combo[i][1] for i in range(3)]
                print(f"       sequence spacing: {spacing}"
                      f"{'   -> CxxC...CxxC zinc-ribbon pattern' if spacing[0] == 3 and spacing[2] == 3 else ''}")
                mean_pl = np.mean([bfac.get((k[0], k[1]), 0.0) for k in combo])
                print(f"       mean pLDDT of the four cysteines: {mean_pl:.1f}")
        if not found:
            print("    none")
            print(f"    (note: a real site can be missed if the model is poor; "
                  f"raise --max-ss to look wider)")

    if args.regions:
        print(f"\n  mean pLDDT (B-factor) by region:")
        for spec in args.regions:
            if ":" in spec:
                rng, label = spec.split(":", 1)
            else:
                rng, label = spec, spec
            try:
                lo, hi = (int(x) for x in rng.split("-"))
            except ValueError:
                print(f"    bad region '{spec}', expected LO-HI[:NAME]", file=sys.stderr)
                continue
            vals = [v for (ch, r), v in bfac.items() if lo <= r <= hi]
            if vals:
                print(f"    {label:<24} {np.mean(vals):5.1f}   "
                      f"(min {min(vals):.1f}, n={len(vals)})")
            else:
                print(f"    {label:<24} -- no residues in range")
    print()


if __name__ == "__main__":
    main()
