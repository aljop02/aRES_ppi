#!/usr/bin/env python3
"""
zn_coordination.py -- measure the coordination geometry of a metal ion in a model.

Companion to zn_site_check.py. That script asks "does this apo model have metal-site
geometry?" by measuring SG...SG distances with nothing in the pocket. This one asks
the complementary question: given a model that *does* contain a metal, is the metal
actually coordinated, and is the geometry chemically correct?

That distinction matters. Handing AlphaFold3 a Zn2+ guarantees a Zn appears in the
output -- it does not guarantee the Zn went anywhere sensible. The diagnostic is
whether the first coordination shell is the expected one and whether bond lengths and
angles match the ideal:

    Zn-S(Cys)        2.30 A
    Zn-N(His NE2)    2.05 A
    S-Zn-S angle   109.47 deg  (tetrahedral)

A site that is real comes out distorted but close. A metal dumped on the surface
comes out with one or two contacts and no angular agreement.

Reads mmCIF (AlphaFold3) and PDB. The B-factor column is pLDDT for AlphaFold output,
including for the ligand, so the per-atom confidence of the metal itself is reported.

Dependencies: numpy.

Usage
-----
  python3 zn_coordination.py model.cif
  python3 zn_coordination.py model.pdb --metal ZN --shell 3.0
"""

from __future__ import annotations

import argparse
import itertools
import math
import sys

import numpy as np

IDEAL = {"SG": 2.30, "NE2": 2.05, "ND1": 2.05, "OD1": 1.95, "OE1": 1.95}
TETRAHEDRAL = math.degrees(math.acos(-1.0 / 3.0))  # 109.47 deg


def load(path: str) -> list[dict]:
    """Parse ATOM/HETATM records from mmCIF or PDB into a flat atom list."""
    atoms, is_cif = [], path.lower().endswith((".cif", ".mmcif"))
    with open(path, errors="replace") as fh:
        for ln in fh:
            if not ln.startswith(("ATOM ", "ATOM  ", "HETATM")):
                if not is_cif and ln.startswith("ENDMDL"):
                    break
                continue
            try:
                if is_cif:
                    # AF3 atom_site order: group id symbol atom alt comp asym entity seq ins x y z occ B
                    f = ln.split()
                    atom, comp, chain = f[3], f[5], f[6]
                    # label_seq_id is "." for non-polymer entities (ligands, ions);
                    # auth_seq_id is the last field. Dropping the row here would
                    # silently discard the metal, which is the whole point.
                    seq = int(f[8]) if f[8].lstrip("-").isdigit() else int(f[-1])
                    xyz = np.array([float(f[10]), float(f[11]), float(f[12])])
                    b = float(f[14])
                else:
                    atom, comp, chain = ln[12:16].strip(), ln[17:20].strip(), ln[21]
                    seq = int(ln[22:26])
                    xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
                    b = float(ln[60:66])
            except (ValueError, IndexError):
                continue
            atoms.append({"atom": atom.upper(), "comp": comp.upper(), "chain": chain,
                          "seq": seq, "xyz": xyz, "b": b})
    return atoms


def main() -> None:
    ap = argparse.ArgumentParser(description="Measure metal coordination geometry.")
    ap.add_argument("structure", help="PDB or mmCIF file")
    ap.add_argument("--metal", default="ZN", help="metal residue name (default ZN)")
    ap.add_argument("--shell", type=float, default=3.0,
                    help="first-shell cutoff in A (default 3.0)")
    args = ap.parse_args()

    atoms = load(args.structure)
    if not atoms:
        sys.exit(f"no atoms parsed from {args.structure}")
    metals = [a for a in atoms if a["comp"] == args.metal.upper()]
    protein = [a for a in atoms if a["comp"] != args.metal.upper()]

    print(f"\n{args.structure}")
    print(f"  atoms {len(atoms)},  {args.metal} ions {len(metals)}")
    if not metals:
        print(f"  no {args.metal} in this model -- use zn_site_check.py for apo models\n")
        return

    for n, m in enumerate(metals, 1):
        print(f"\n  {args.metal} #{n}   pLDDT {m['b']:.1f}")
        dists = sorted(((float(np.linalg.norm(m["xyz"] - a["xyz"])), a) for a in protein),
                       key=lambda t: t[0])
        shell = [(d, a) for d, a in dists if d <= args.shell]

        print(f"    first coordination shell (< {args.shell} A): "
              f"{len(shell)} atom{'s' if len(shell) != 1 else ''}")
        if not shell:
            print("    -> UNCOORDINATED. The metal is not in a site.")
            print(f"       nearest atom: {dists[0][1]['comp']}{dists[0][1]['seq']} "
                  f"{dists[0][1]['atom']} at {dists[0][0]:.2f} A")
            continue

        for d, a in shell:
            ideal = IDEAL.get(a["atom"])
            dev = f"  (ideal {ideal:.2f}, {d - ideal:+.2f})" if ideal else ""
            print(f"      {d:5.2f} A   {a['comp']}{a['seq']:<4} {a['atom']:<4} "
                  f"pLDDT {a['b']:5.1f}{dev}")

        bonded = [a for _, a in shell if a["atom"] in IDEAL]
        if len(bonded) >= 2:
            print(f"\n    ligand-{args.metal}-ligand angles "
                  f"(ideal tetrahedral {TETRAHEDRAL:.2f} deg):")
            angles = []
            for a, b in itertools.combinations(bonded, 2):
                v1, v2 = a["xyz"] - m["xyz"], b["xyz"] - m["xyz"]
                cos = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
                ang = math.degrees(math.acos(max(-1.0, min(1.0, cos))))
                angles.append(ang)
                print(f"      {a['comp']}{a['seq']}-{args.metal}-{b['comp']}{b['seq']}: "
                      f"{ang:6.1f}")
            if len(bonded) == 4:
                mean = float(np.mean(angles))
                print(f"      mean {mean:.1f} deg  "
                      f"(deviation from tetrahedral {abs(mean - TETRAHEDRAL):.1f} deg, "
                      f"spread {max(angles) - min(angles):.1f} deg)")

        # Is anything else nearby that could have competed for the metal?
        others = [(d, a) for d, a in dists
                  if d > args.shell and a["atom"] in IDEAL][:3]
        if others:
            print(f"\n    next candidate ligands not in the shell:")
            for d, a in others:
                print(f"      {d:6.2f} A   {a['comp']}{a['seq']} {a['atom']}")
    print()


if __name__ == "__main__":
    main()
