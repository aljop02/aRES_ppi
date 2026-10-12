#!/usr/bin/env python3
"""
make_haddock_restraints.py -- build HADDOCK ambiguous interaction restraints (AIRs).

Writes a .tbl file for HADDOCK 2.4 / 3 from a residue selection, filtering the
selection down to solvent-exposed residues first (buried residues cannot be part
of an interface, and restraining them drags the docking somewhere wrong).

The design this project needs is the ONE-SIDED run: restrain only the epitope you
know experimentally (gp5 residues ~100-125, defined by the escape mutations), and
leave the aRES side completely free. Whichever aRES domain HADDOCK then picks is a
genuine prediction rather than an assumption you fed in.

Relative SASA uses Shrake-Rupley (Biopython) over the ISOLATED chain, divided by the
Tien et al. 2013 theoretical maximum ASA per residue type. Compute it on the unbound
monomer, not on a complex -- in a complex the interface is already buried and
everything you care about reads as inaccessible.

Dependencies: biopython, numpy.

Usage
-----
  # One-sided: gp5 loop restrained, aRES free. THE run to start.
  python3 make_haddock_restraints.py \
      --mol-a gp5.pdb      --seg-a A --active-a 100-125 \
      --mol-b ares.pdb     --seg-b B --passive-b surface \
      --out air_onesided.tbl

  # Two-sided: both interfaces asserted (use only to refine a pose you trust)
  python3 make_haddock_restraints.py \
      --mol-a gp5.pdb  --seg-a A --active-a 100-125 \
      --mol-b ares.pdb --seg-b B --active-b 95-115 \
      --out air_twosided.tbl

  # Just report accessibility, write nothing
  python3 make_haddock_restraints.py --mol-a ares.pdb --seg-a A --report 1-120
"""

from __future__ import annotations

import argparse
import sys

from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.SASA import ShrakeRupley

# Tien et al. 2013, theoretical maximum accessible surface area (A^2) per residue.
# doi:10.1371/journal.pone.0080635
MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLU": 223.0, "GLN": 225.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}


def parse_selection(spec: str) -> list[int]:
    """'100-125,130,140-145' -> [100..125, 130, 140..145]"""
    out: list[int] = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part.lstrip("-"):
            lo_s, hi_s = part.split("-", 1)
            try:
                lo, hi = int(lo_s), int(hi_s)
            except ValueError:
                sys.exit(f"bad range '{part}'")
            if lo > hi:
                sys.exit(f"range '{part}' is reversed")
            out.extend(range(lo, hi + 1))
        else:
            try:
                out.append(int(part))
            except ValueError:
                sys.exit(f"bad residue '{part}'")
    return sorted(set(out))


def load_chain(path: str, chain_id: str | None):
    parser = MMCIFParser(QUIET=True) if path.lower().endswith((".cif", ".mmcif")) \
        else PDBParser(QUIET=True)
    structure = parser.get_structure("m", path)
    model = next(structure.get_models())
    chains = list(model)
    if chain_id:
        for ch in chains:
            if ch.id == chain_id:
                return ch
        sys.exit(f"{path}: chain '{chain_id}' not found (have: "
                 f"{[c.id for c in chains]})")
    if len(chains) > 1:
        print(f"  note: {path} has chains {[c.id for c in chains]}; "
              f"using '{chains[0].id}'. Pass --chain-a/--chain-b to choose.",
              file=sys.stderr)
    return chains[0]


def relative_sasa(chain) -> dict[int, tuple[str, float]]:
    """{resseq: (resname, relative SASA in percent)} for the isolated chain."""
    # Strip waters/hetero so SASA is computed on the polymer alone.
    for res in [r for r in chain if r.id[0] != " "]:
        chain.detach_child(res.id)
    ShrakeRupley().compute(chain, level="R")
    out = {}
    for res in chain:
        name = res.get_resname().upper()
        ref = MAX_ASA.get(name)
        if ref is None:
            continue  # non-standard residue; no reference value
        out[res.id[1]] = (name, 100.0 * res.sasa / ref)
    return out


def filter_exposed(selection, sasa, threshold, label):
    keep, dropped, missing = [], [], []
    for r in selection:
        if r not in sasa:
            missing.append(r)
        elif sasa[r][1] >= threshold:
            keep.append(r)
        else:
            dropped.append((r, sasa[r][0], sasa[r][1]))
    if missing:
        print(f"  ! {label}: {len(missing)} selected residue(s) not in the structure: "
              f"{', '.join(map(str, missing[:15]))}"
              f"{' ...' if len(missing) > 15 else ''}", file=sys.stderr)
    if dropped:
        print(f"  {label}: dropped {len(dropped)} buried residue(s) "
              f"(< {threshold}% rel. SASA):")
        for r, nm, v in dropped:
            print(f"      {nm}{r}  {v:.1f}%")
    return keep


def write_air(fh, actives_a, seg_a, targets_b, seg_b, lower=2.0, upper=0.0, dist=2.0):
    """One block per active residue, restrained ambiguously to every target."""
    n = 0
    for a in actives_a:
        fh.write(f"assign ( resid {a:>4d} and segid {seg_a} )\n")
        fh.write("       (\n")
        for i, b in enumerate(targets_b):
            if i:
                fh.write("     or\n")
            fh.write(f"        ( resid {b:>4d} and segid {seg_b} )\n")
        fh.write(f"       )  {dist} {lower} {upper}\n\n")
        n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Build HADDOCK AIR restraints, filtered by solvent accessibility.")
    ap.add_argument("--mol-a", required=True, help="PDB/mmCIF of molecule A")
    ap.add_argument("--chain-a", help="chain to read from --mol-a")
    ap.add_argument("--seg-a", default="A", help="HADDOCK segid for A (default A)")
    ap.add_argument("--active-a", help="active residues on A, e.g. 100-125")
    ap.add_argument("--passive-a", help="passive residues on A, or 'surface'")

    ap.add_argument("--mol-b", help="PDB/mmCIF of molecule B")
    ap.add_argument("--chain-b", help="chain to read from --mol-b")
    ap.add_argument("--seg-b", default="B", help="HADDOCK segid for B (default B)")
    ap.add_argument("--active-b", help="active residues on B")
    ap.add_argument("--passive-b", help="passive residues on B, or 'surface'")

    ap.add_argument("--min-sasa", type=float, default=40.0,
                    help="relative SASA %% threshold for an active residue (default 40)")
    ap.add_argument("--surface-sasa", type=float, default=15.0,
                    help="relative SASA %% threshold for 'surface' passives (default 15)")
    ap.add_argument("--report", help="just print accessibility for these residues of A")
    ap.add_argument("--out", help="output .tbl path")
    args = ap.parse_args()

    chain_a = load_chain(args.mol_a, args.chain_a)
    sasa_a = relative_sasa(chain_a)
    print(f"molecule A: {args.mol_a} chain {chain_a.id}, {len(sasa_a)} standard residues")

    if args.report:
        sel = parse_selection(args.report)
        print(f"\n  {'residue':>12}   rel.SASA   exposed(>={args.min_sasa}%)")
        for r in sel:
            if r not in sasa_a:
                print(f"  {r:>12}   --  not in structure")
                continue
            nm, v = sasa_a[r]
            print(f"  {nm + str(r):>12}   {v:6.1f}%   {'yes' if v >= args.min_sasa else 'no'}")
        return

    if not args.mol_b:
        sys.exit("--mol-b is required unless --report is used")
    if not args.out:
        sys.exit("--out is required")
    if not (args.active_a or args.active_b):
        sys.exit("at least one of --active-a / --active-b is required")

    chain_b = load_chain(args.mol_b, args.chain_b)
    sasa_b = relative_sasa(chain_b)
    print(f"molecule B: {args.mol_b} chain {chain_b.id}, {len(sasa_b)} standard residues")

    def resolve(spec, sasa, label, *, as_active):
        if not spec:
            return []
        if spec.lower() == "surface":
            sel = sorted(r for r, (_, v) in sasa.items() if v >= args.surface_sasa)
            print(f"  {label}: 'surface' -> {len(sel)} residues "
                  f"(>= {args.surface_sasa}% rel. SASA)")
            return sel
        sel = parse_selection(spec)
        thr = args.min_sasa if as_active else args.surface_sasa
        return filter_exposed(sel, sasa, thr, label)

    print()
    act_a = resolve(args.active_a, sasa_a, "active A", as_active=True)
    pas_a = resolve(args.passive_a, sasa_a, "passive A", as_active=False)
    act_b = resolve(args.active_b, sasa_b, "active B", as_active=True)
    pas_b = resolve(args.passive_b, sasa_b, "passive B", as_active=False)

    if args.active_a and not act_a:
        sys.exit("every selected active residue on A is buried -- check the selection "
                 "and that the SASA was computed on the UNBOUND monomer")
    if args.active_b and not act_b:
        sys.exit("every selected active residue on B is buried -- same check")

    targets_b = sorted(set(act_b) | set(pas_b))
    targets_a = sorted(set(act_a) | set(pas_a))
    if act_a and not targets_b:
        sys.exit("A has active residues but B has no active or passive targets; "
                 "pass --passive-b surface for a one-sided run")
    if act_b and not targets_a:
        sys.exit("B has active residues but A has no active or passive targets")

    with open(args.out, "w") as fh:
        fh.write("! HADDOCK ambiguous interaction restraints\n")
        fh.write(f"! generated by make_haddock_restraints.py\n")
        fh.write(f"! active A ({args.seg_a}): {args.active_a or '-'}\n")
        fh.write(f"! active B ({args.seg_b}): {args.active_b or '-'}\n")
        fh.write(f"! min rel. SASA for actives: {args.min_sasa}%\n!\n")
        n = 0
        if act_a:
            n += write_air(fh, act_a, args.seg_a, targets_b, args.seg_b)
        if act_b:
            n += write_air(fh, act_b, args.seg_b, targets_a, args.seg_a)

    sided = "two-sided" if (act_a and act_b) else "ONE-SIDED"
    print(f"\nwrote {args.out}")
    print(f"  {n} AIR blocks, {sided}")
    print(f"  A: {len(act_a)} active, {len(pas_a)} passive")
    print(f"  B: {len(act_b)} active, {len(pas_b)} passive")
    if sided == "ONE-SIDED":
        free = args.seg_b if act_a else args.seg_a
        print(f"  -> segid {free} is unrestrained: whatever interface it forms is a "
              f"prediction, not an assumption.")


if __name__ == "__main__":
    main()
