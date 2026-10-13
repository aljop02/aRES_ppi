#!/usr/bin/env python3
"""
prepare_haddock_input.py -- make a PDB safe to hand to the HADDOCK 2.4 web server.

Two jobs, both of which fix real failure modes seen in this project:

1. **Stamp an explicit segid** (columns 73-76). AlphaFold models carry a chain ID
   but a blank segid, and every restraint we write selects on `segid`. Letting the
   server infer the segid from the chain is an unnecessary dependency -- and an
   empty restraint selection is exactly what silently killed job 780457.

2. **Insert a structural Zn(2+) into a Cys4 site** and emit the unambiguous
   restraints that hold it there. aRES has a Zn(Cys)4 ribbon at C8/C11/C31/C34
   (docs/07) that was absent from every docking run so far, leaving four
   uncoordinated thiolates in a tight cluster for flexible refinement to distort.

HADDOCK ion conventions (residue `ZN2`, atom `ZN+2`, ion sharing the protein's
chain so it can move with it during rigid-body docking, and the atom name quoted
in restraint files because of the `+`) come from the BioExcel/HADDOCK user forum,
not from a spec this script can check. VERIFY the naming against the server's own
supported-ions page before spending queue time -- see docs/10 step 2.

The zinc is placed at the point that best equidistributes the Zn-SG distances
around the target bond length, which for a near-ideal tetrahedron is essentially
the SG centroid.

Usage
-----
  # segid only
  python3 prepare_haddock_input.py in.pdb -o out.pdb --segid B

  # segid + zinc site, writing the coordination restraints alongside
  python3 prepare_haddock_input.py aRES.pdb -o aRES_Zn.pdb --segid A \
      --zinc-cys 8,11,31,34 --zinc-resid 500 --tbl zn_unambig.tbl
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

ZN_SG_TARGET = 2.30   # Angstrom, Zn(2+)-S(cys) bond
ZN_SG_TOL = 0.10      # +/- on the coordination restraint
SG_SG_TOL = 0.20      # +/- on the tetrahedron edge restraints


def read_pdb(path: str) -> list[str]:
    with open(path, encoding="utf-8") as fh:
        return fh.read().splitlines()


def pad(line: str, n: int = 80) -> str:
    return line.ljust(n)


def stamp_segid(lines: list[str], segid: str) -> list[str]:
    """Write segid into columns 73-76 of every ATOM/HETATM record."""
    seg = segid.ljust(4)[:4]
    out = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            line = pad(line)
            line = line[:72] + seg + line[76:]
        out.append(line.rstrip())
    return out


def sg_coords(lines: list[str], resids: list[int], chain: str | None):
    """Return {resid: xyz} for the SG atom of each named cysteine."""
    found: dict[int, np.ndarray] = {}
    for line in lines:
        if not line.startswith("ATOM"):
            continue
        line = pad(line)
        if line[12:16].strip() != "SG":
            continue
        if chain and line[21] != chain:
            continue
        try:
            r = int(line[22:26])
        except ValueError:
            continue
        if r in resids:
            found[r] = np.array([float(line[30:38]), float(line[38:46]),
                                 float(line[46:54])])
    missing = [r for r in resids if r not in found]
    if missing:
        raise SystemExit(f"error: no SG atom found for residue(s) {missing}")
    return found


def place_zinc(sgs: dict[int, np.ndarray]) -> np.ndarray:
    """Point minimising sum((|Zn-SG_i| - target)^2); start from the centroid."""
    P = np.array(list(sgs.values()))
    z = P.mean(axis=0)
    for _ in range(200):                       # a few Gauss-Newton style steps
        d = np.linalg.norm(P - z, axis=1)
        if np.all(d < 1e-9):
            break
        grad = ((d - ZN_SG_TARGET)[:, None] * (z - P) / d[:, None]).sum(axis=0)
        z = z - 0.5 * grad / len(P)
    return z


def rename_cys(lines: list[str], resids: list[int], chain: str | None,
               new: str = "CYF") -> tuple[list[str], int]:
    """CYS -> CYF for the coordinating cysteines (thiolate, no SH proton)."""
    out, n = [], 0
    for line in lines:
        if line.startswith("ATOM"):
            line = pad(line)
            if line[17:20].strip() == "CYS" and (not chain or line[21] == chain):
                try:
                    r = int(line[22:26])
                except ValueError:
                    r = None
                if r in resids:
                    line = line[:17] + new.ljust(3) + line[20:]
                    n += 1
        out.append(line.rstrip())
    return out, n


def zinc_record(serial: int, resid: int, chain: str, segid: str,
                xyz: np.ndarray) -> str:
    """HETATM for HADDOCK: 4-char atom name ZN+2 starts in column 13."""
    rec = (
        f"HETATM{serial:>5d} "          # 1-6 record, 7-11 serial, 12 space
        f"{'ZN+2':<4s}"                  # 13-16 atom name (4 chars, no pad left)
        f" "                             # 17 altLoc
        f"{'ZN2':>3s}"                   # 18-20 resName
        f" "                             # 21 space
        f"{chain:1s}"                    # 22 chainID
        f"{resid:>4d}"                   # 23-26 resSeq
        f"    "                          # 27 iCode + 28-30 pad
        f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
        f"{1.00:6.2f}{0.00:6.2f}"        # occupancy, B
        f"      "                        # 67-72
        f"{segid:<4s}"                   # 73-76 segid
        f"{'ZN':>2s}"                    # 77-78 element
    )
    return rec


def write_tbl(path: str, segid: str, zinc_resid: int,
              sgs: dict[int, np.ndarray], zpos: np.ndarray) -> int:
    lines = [
        "! Unambiguous restraints: structural Zn(Cys)4 site",
        "! Generated by analysis/scripts/prepare_haddock_input.py",
        f"! Zn is segid {segid} resid {zinc_resid}, atom \"ZN+2\"",
        "! 4 Zn-SG coordination restraints + 6 SG-SG tetrahedron edges",
        "!",
    ]
    n = 0
    for r, xyz in sorted(sgs.items()):
        d = float(np.linalg.norm(xyz - zpos))
        lines.append(f"! Zn-SG {r}: modelled {d:.2f} A")
        lines.append(f"assign ( segid {segid} and resid {r:>4d} and name SG )")
        lines.append(f"       ( segid {segid} and resid {zinc_resid:>4d} "
                     f"and name \"ZN+2\" )  "
                     f"{ZN_SG_TARGET:.2f} {ZN_SG_TOL:.2f} {ZN_SG_TOL:.2f}")
        n += 1
    rs = sorted(sgs)
    for i in range(len(rs)):
        for j in range(i + 1, len(rs)):
            a, b = rs[i], rs[j]
            d = float(np.linalg.norm(sgs[a] - sgs[b]))
            lines.append(f"! SG{a}-SG{b} tetrahedron edge, modelled {d:.2f} A")
            lines.append(f"assign ( segid {segid} and resid {a:>4d} and name SG )")
            lines.append(f"       ( segid {segid} and resid {b:>4d} and name SG )  "
                         f"{d:.2f} {SG_SG_TOL:.2f} {SG_SG_TOL:.2f}")
            n += 1
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdb", help="input PDB")
    ap.add_argument("-o", "--out", required=True, help="output PDB")
    ap.add_argument("--segid", required=True, help="segid to stamp (e.g. A)")
    ap.add_argument("--chain", help="only touch this chain ID (default: all)")
    ap.add_argument("--zinc-cys", help="comma-separated cysteine resids, e.g. 8,11,31,34")
    ap.add_argument("--zinc-resid", type=int, default=500, help="resid for the Zn (default 500)")
    ap.add_argument("--tbl", help="write the Zn restraints here")
    ap.add_argument("--no-cyf", action="store_true",
                    help="do not rename coordinating CYS to CYF")
    args = ap.parse_args()

    lines = read_pdb(args.pdb)
    chain = args.chain

    if args.zinc_cys:
        resids = [int(x) for x in args.zinc_cys.split(",") if x.strip()]
        if len(resids) != 4:
            print(f"warning: expected 4 cysteines, got {len(resids)}", file=sys.stderr)
        sgs = sg_coords(lines, resids, chain)
        zpos = place_zinc(sgs)

        print("Zn(Cys)4 site:", file=sys.stderr)
        for r, xyz in sorted(sgs.items()):
            print(f"  Zn-SG {r:>4d}: {np.linalg.norm(xyz - zpos):.3f} A", file=sys.stderr)
        edges = [float(np.linalg.norm(sgs[a] - sgs[b]))
                 for i, a in enumerate(sorted(sgs)) for b in sorted(sgs)[i + 1:]]
        print(f"  SG-SG edges: {min(edges):.2f}-{max(edges):.2f} A "
              f"(mean {sum(edges)/len(edges):.2f})", file=sys.stderr)
        print(f"  Zn placed at {zpos[0]:.3f} {zpos[1]:.3f} {zpos[2]:.3f}", file=sys.stderr)

        if not args.no_cyf:
            lines, nren = rename_cys(lines, resids, chain)
            print(f"  renamed {nren} atoms CYS -> CYF", file=sys.stderr)

        # find the chain to attach the ion to, and the highest serial
        zn_chain = chain
        serial = 0
        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                p = pad(line)
                try:
                    serial = max(serial, int(p[6:11]))
                except ValueError:
                    pass
                if zn_chain is None:
                    zn_chain = p[21]
        used = {int(pad(l)[22:26]) for l in lines
                if l.startswith(("ATOM", "HETATM")) and pad(l)[21] == zn_chain
                and pad(l)[22:26].strip().lstrip("-").isdigit()}
        if args.zinc_resid in used:
            raise SystemExit(f"error: resid {args.zinc_resid} already used in chain "
                             f"{zn_chain}; pick another with --zinc-resid")

        rec = zinc_record(serial + 1, args.zinc_resid, zn_chain, args.segid, zpos)
        # insert before any trailing TER/END
        body = [l for l in lines if not l.startswith(("TER", "END"))]
        body.append(rec)
        body.append("TER")
        body.append("END")
        lines = body

        if args.tbl:
            n = write_tbl(args.tbl, args.segid, args.zinc_resid, sgs, zpos)
            print(f"  wrote {args.tbl}: {n} restraints", file=sys.stderr)

    lines = stamp_segid(lines, args.segid)
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    natom = sum(1 for l in lines if l.startswith("ATOM"))
    nhet = sum(1 for l in lines if l.startswith("HETATM"))
    print(f"wrote {args.out}: {natom} ATOM, {nhet} HETATM, segid {args.segid!r}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
