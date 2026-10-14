#!/usr/bin/env python3
"""
homolog_compare.py -- map residue numbering and superpose across aRES orthologs.

The problem this solves: every residue number in this project (C8/C11/C31/C34,
E89, E93, D104) is in the *UC5* frame, and orthologs have different lengths --
Klebsiella aRES is 443 aa against UC5's 459. Quoting "D104" for a Klebsiella model
is meaningless until the frames are aligned, and eyeballing an alignment is how
frame errors get into documents (see docs/06 on the UC5/S_34AVR mix-up).

So: global pairwise alignment (Needleman-Wunsch with affine gaps, BLOSUM62) between
the reference and each ortholog, then

  * per-region percent identity
  * what each reference position of interest maps to in every ortholog
  * CA superposition (Kabsch) over a chosen reference region, using only
    alignment-equivalent pairs, with RMSD

Sequences are read from the structures themselves, not from a separate FASTA, so
the alignment and the coordinates cannot drift apart.

Caveat on the RMSD: it is computed over aligned pairs after a least-squares fit, so
it answers "is this the same fold in the same conformation" and nothing more. A low
RMSD over a region says the models agree there; it does not make either one right.

Dependencies: numpy.

Usage
-----
  python3 homolog_compare.py ref.cif other1.cif other2.cif \
      --names UC5 E.chengduensis Klebsiella \
      --regions 1-48:Zn-ribbon 49-120:NTD1 126-447:RES \
      --map 8 11 31 34 89 93 104 \
      --super 1-120
"""

from __future__ import annotations

import argparse
import itertools
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from af3_zn_site import parse_cif_atoms  # noqa: E402  (same directory)

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

AA = "ARNDCQEGHILKMFPSTWYV"

_B62 = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""

IDX = {c: i for i, c in enumerate(AA)}
BLOSUM62 = np.array([[int(x) for x in row.split()] for row in _B62.strip().split("\n")])

GAP_OPEN = -11
GAP_EXTEND = -1


def score(a: str, b: str) -> int:
    ia, ib = IDX.get(a), IDX.get(b)
    if ia is None or ib is None:
        return -1
    return int(BLOSUM62[ia, ib])


# ------------------------------------------------------------------ alignment


def align(s1: str, s2: str):
    """Needleman-Wunsch with affine gaps (Gotoh). -> list of (i, j) 0-based pairs.

    Gaps are represented by None on one side.
    """
    n, m = len(s1), len(s2)
    NEG = -10**9
    # M: match state, X: gap in s2 (consume s1), Y: gap in s1 (consume s2)
    M = np.full((n + 1, m + 1), NEG, dtype=np.int64)
    X = np.full((n + 1, m + 1), NEG, dtype=np.int64)
    Y = np.full((n + 1, m + 1), NEG, dtype=np.int64)
    pM = np.zeros((n + 1, m + 1), dtype=np.int8)
    pX = np.zeros((n + 1, m + 1), dtype=np.int8)
    pY = np.zeros((n + 1, m + 1), dtype=np.int8)

    M[0, 0] = 0
    for i in range(1, n + 1):
        X[i, 0] = GAP_OPEN + (i - 1) * GAP_EXTEND
        pX[i, 0] = 1
    for j in range(1, m + 1):
        Y[0, j] = GAP_OPEN + (j - 1) * GAP_EXTEND
        pY[0, j] = 2

    for i in range(1, n + 1):
        a = s1[i - 1]
        srow = np.array([score(a, b) for b in s2], dtype=np.int64)
        for j in range(1, m + 1):
            s = srow[j - 1]
            cand = (M[i - 1, j - 1], X[i - 1, j - 1], Y[i - 1, j - 1])
            k = int(np.argmax(cand))
            M[i, j] = cand[k] + s
            pM[i, j] = k

            cand = (M[i - 1, j] + GAP_OPEN, X[i - 1, j] + GAP_EXTEND)
            k = int(np.argmax(cand))
            X[i, j] = cand[k]
            pX[i, j] = 0 if k == 0 else 1

            cand = (M[i, j - 1] + GAP_OPEN, Y[i, j - 1] + GAP_EXTEND)
            k = int(np.argmax(cand))
            Y[i, j] = cand[k]
            pY[i, j] = 0 if k == 0 else 2

    i, j = n, m
    state = int(np.argmax((M[n, m], X[n, m], Y[n, m])))
    pairs = []
    while i > 0 or j > 0:
        if state == 0:
            pairs.append((i - 1, j - 1))
            state = int(pM[i, j])
            i, j = i - 1, j - 1
        elif state == 1:
            pairs.append((i - 1, None))
            state = int(pX[i, j])
            i -= 1
        else:
            pairs.append((None, j - 1))
            state = int(pY[i, j])
            j -= 1
        if i == 0 and j == 0:
            break
        if i == 0 and state != 2:
            state = 2
        if j == 0 and state != 1:
            state = 1
    return list(reversed(pairs))


# ----------------------------------------------------------------- structures


def read_chain(path: str):
    """-> (seq, resnums, ca_coords, plddt_by_resnum) for the longest protein chain."""
    atoms = parse_cif_atoms(path)
    by_chain = {}
    for a in atoms:
        if a["resname"] in THREE_TO_ONE:
            by_chain.setdefault(a["chain"], []).append(a)
    if not by_chain:
        raise SystemExit(f"{path}: no protein residues")
    chain = max(by_chain, key=lambda c: len({a["resseq"] for a in by_chain[c]}))
    sel = by_chain[chain]

    resnames, cas, plddt = {}, {}, {}
    for a in sel:
        resnames.setdefault(a["resseq"], a["resname"])
        plddt.setdefault(a["resseq"], []).append(a["plddt"])
        if a["atom"] == "CA":
            cas[a["resseq"]] = a["xyz"]
    nums = sorted(resnames)
    seq = "".join(THREE_TO_ONE[resnames[r]] for r in nums)
    ca = np.array([cas.get(r, np.full(3, np.nan)) for r in nums])
    pl = {r: float(np.nanmean(v)) for r, v in plddt.items()}
    return seq, nums, ca, pl


def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """RMSD after optimal superposition of P onto Q."""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt
    Pr = Pc @ R
    return float(np.sqrt(np.mean(np.sum((Pr - Qc) ** 2, axis=1))))


def parse_regions(specs):
    out = []
    for spec in specs:
        rng, _, label = spec.partition(":")
        lo, hi = (int(x) for x in rng.split("-"))
        out.append((label or rng, lo, hi))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Compare aRES orthologs by sequence and structure.")
    ap.add_argument("structures", nargs="+", help="first is the reference frame")
    ap.add_argument("--names", nargs="*", default=[])
    ap.add_argument("--regions", nargs="*", default=[], metavar="LO-HI:NAME")
    ap.add_argument("--map", nargs="*", type=int, default=[], metavar="RESNUM",
                    help="reference residue numbers to map onto each ortholog")
    ap.add_argument("--super", nargs="*", default=[], metavar="LO-HI:NAME",
                    help="reference regions to superpose over (also whole chain)")
    ap.add_argument("--profile", help="write a windowed identity profile TSV here")
    ap.add_argument("--window", type=int, default=15,
                    help="window size for the identity profile (default 15)")
    ap.add_argument("--conserved", nargs="*", default=[], metavar="LO-HI:NAME",
                    help="regions to list fully invariant positions for")
    ap.add_argument("--min-indel", type=int, default=2,
                    help="shortest indel run to report (default 2)")
    args = ap.parse_args()

    names = list(args.names) + [
        os.path.basename(p) for p in args.structures[len(args.names):]
    ]
    chains = [read_chain(p) for p in args.structures]
    ref_name, (ref_seq, ref_nums, ref_ca, ref_pl) = names[0], chains[0]

    print("=" * 78)
    print(f"reference frame: {ref_name}  ({len(ref_seq)} aa)")
    print("=" * 78)
    for nm, (seq, _, _, _) in zip(names, chains):
        print(f"  {nm:<20} {len(seq)} aa")

    regions = parse_regions(args.regions)
    super_regions = parse_regions(args.super)

    alignments = {}
    for nm, (seq, nums, ca, pl) in list(zip(names, chains))[1:]:
        pairs = align(ref_seq, seq)
        alignments[nm] = (pairs, seq, nums, ca, pl)

    # ---- identity, overall and per region (in reference numbering)
    print(f"\npercent identity to {ref_name}")
    hdr = f"  {'ortholog':<20}{'overall':>9}" + "".join(f"{lab:>14}" for lab, _, _ in regions)
    print(hdr)
    for nm, (pairs, seq, nums, ca, pl) in alignments.items():
        def ident(lo=None, hi=None):
            n = ok = 0
            for i, j in pairs:
                if i is None or j is None:
                    continue
                rn = ref_nums[i]
                if lo is not None and not (lo <= rn <= hi):
                    continue
                n += 1
                ok += ref_seq[i] == seq[j]
            return (100.0 * ok / n) if n else float("nan")

        line = f"  {nm:<20}{ident():>8.1f}%"
        for _, lo, hi in regions:
            line += f"{ident(lo, hi):>13.1f}%"
        print(line)

    # ---- map positions of interest
    if args.map:
        print(f"\nreference positions mapped onto each ortholog")
        print(f"  {'ref':<10}" + "".join(f"{nm:>22}" for nm in list(alignments)))
        for rn in args.map:
            try:
                i = ref_nums.index(rn)
            except ValueError:
                print(f"  {rn} not in reference")
                continue
            cells = []
            for nm, (pairs, seq, nums, ca, pl) in alignments.items():
                j = next((jj for ii, jj in pairs if ii == i), None)
                if j is None:
                    cells.append("-- gap --")
                else:
                    same = "=" if seq[j] == ref_seq[i] else " "
                    cells.append(f"{seq[j]}{nums[j]}{same} pLDDT {pl.get(nums[j], float('nan')):.0f}")
            print(f"  {ref_seq[i]}{rn:<9}" + "".join(f"{c:>22}" for c in cells))
        print("  ('=' marks an identical residue; pLDDT is the ortholog's own)")

    # ---- positions invariant across every ortholog, within a region
    # In a module that is otherwise diverging, the positions that do not move are
    # the ones under constraint, and so the only defensible mutagenesis targets.
    for lab, lo, hi in parse_regions(args.conserved):
        rows = []
        for i, rn in enumerate(ref_nums):
            if not (lo <= rn <= hi):
                continue
            aas, ok = [], True
            for nm, (pairs, seq, nums, ca, pl) in alignments.items():
                j = next((jj for ii, jj in pairs if ii == i), None)
                if j is None or seq[j] != ref_seq[i]:
                    ok = False
                    break
                aas.append(f"{nums[j]}")
            if ok:
                rows.append((ref_seq[i], rn, aas))
        n_region = sum(1 for rn in ref_nums if lo <= rn <= hi)
        print(f"\ninvariant across all {len(alignments) + 1} orthologs in {lab} "
              f"({len(rows)}/{n_region} positions = "
              f"{100.0 * len(rows) / n_region if n_region else float('nan'):.0f}%)")
        print("  " + ", ".join(f"{aa}{rn}" for aa, rn, _ in rows))

    # ---- indels, in reference numbering
    print(f"\ninsertions / deletions relative to {ref_name} "
          f"(runs of {args.min_indel}+ positions)")
    for nm, (pairs, seq, nums, ca, pl) in alignments.items():
        runs, cur, last_ref = [], None, ref_nums[0]
        for i, j in pairs:
            # "del" = present in the reference, absent here; "ins" = the reverse.
            kind = "del" if j is None else ("ins" if i is None else None)
            if kind is None:
                if cur:
                    runs.append(cur)
                    cur = None
                last_ref = ref_nums[i]
                continue
            if cur and cur[0] == kind:
                cur[2] += 1
            else:
                if cur:
                    runs.append(cur)
                cur = [kind, last_ref, 1]
        if cur:
            runs.append(cur)
        runs = [r for r in runs if r[2] >= args.min_indel]
        if not runs:
            print(f"  {nm:<20} none of length {args.min_indel}+")
            continue
        label = []
        for kind, anchor, n in runs:
            which = next((lab for lab, lo, hi in regions if lo <= anchor <= hi), "-")
            label.append(f"{kind} {n} aa at ref~{anchor} [{which}]")
        print(f"  {nm:<20} " + "; ".join(label))

    # ---- superposition
    targets = super_regions or []
    print(f"\nCA RMSD to {ref_name} after Kabsch superposition (aligned pairs only)")
    print(f"  {'ortholog':<20}{'region':<16}{'n_pairs':>9}{'RMSD':>9}")
    for nm, (pairs, seq, nums, ca, pl) in alignments.items():
        specs = [("whole chain", None, None)] + targets
        for lab, lo, hi in specs:
            P, Q = [], []
            for i, j in pairs:
                if i is None or j is None:
                    continue
                rn = ref_nums[i]
                if lo is not None and not (lo <= rn <= hi):
                    continue
                if np.isnan(ref_ca[i]).any() or np.isnan(ca[j]).any():
                    continue
                Q.append(ref_ca[i])
                P.append(ca[j])
            if len(P) < 4:
                print(f"  {nm:<20}{lab:<16}{len(P):>9}{'--':>9}")
                continue
            r = kabsch_rmsd(np.array(P), np.array(Q))
            print(f"  {nm:<20}{lab:<16}{len(P):>9}{r:>9.2f}")

    # ---- inter-domain arrangement
    # If every domain superposes well on its own but the whole chain does not, the
    # modules are conserved and their relative arrangement is not. Measuring that:
    # fit on domain A, then report domain B's RMSD in A's frame. A large number
    # means a different hinge angle, which is a structural claim worth making
    # separately from "the fold is conserved".
    if len(targets) >= 2:
        print(f"\ninter-domain arrangement: fit on one region, measure another in that frame")
        print(f"  {'ortholog':<18}{'fit on':<14}{'measured':<14}{'n':>5}{'RMSD':>9}")
        for nm, (pairs, seq, nums, ca, pl) in alignments.items():
            eq = [(i, j) for i, j in pairs if i is not None and j is not None]
            for fit_lab, fit_lo, fit_hi in targets:
                sel_fit = [(i, j) for i, j in eq if fit_lo <= ref_nums[i] <= fit_hi]
                if len(sel_fit) < 4:
                    continue
                P = np.array([ca[j] for _, j in sel_fit])
                Q = np.array([ref_ca[i] for i, _ in sel_fit])
                Pm, Qm = P.mean(0), Q.mean(0)
                V, S, Wt = np.linalg.svd((P - Pm).T @ (Q - Qm))
                d = np.sign(np.linalg.det(V @ Wt))
                R = V @ np.diag([1.0, 1.0, d]) @ Wt
                for lab2, lo2, hi2 in targets:
                    if lab2 == fit_lab:
                        continue
                    sel = [(i, j) for i, j in eq if lo2 <= ref_nums[i] <= hi2]
                    if len(sel) < 4:
                        continue
                    Pb = np.array([ca[j] for _, j in sel])
                    Qb = np.array([ref_ca[i] for i, _ in sel])
                    Pbr = (Pb - Pm) @ R + Qm
                    rms = float(np.sqrt(np.mean(np.sum((Pbr - Qb) ** 2, axis=1))))
                    print(f"  {nm:<18}{fit_lab:<14}{lab2:<14}{len(sel):>5}{rms:>9.2f}")

    # ---- windowed identity profile
    if args.profile:
        w = args.window
        print(f"\nwindowed identity profile written to {args.profile} (window {w} aa)")
        cols = ["ref_resnum", "ref_aa"] + [f"id_{nm}" for nm in alignments]
        prof = {}
        for nm, (pairs, seq, nums, ca, pl) in alignments.items():
            eq = {i: j for i, j in pairs if i is not None and j is not None}
            per_pos = []
            for i in range(len(ref_seq)):
                j = eq.get(i)
                per_pos.append(1.0 if (j is not None and seq[j] == ref_seq[i]) else 0.0)
            arr = np.array(per_pos)
            half = w // 2
            sm = np.array([
                100.0 * arr[max(0, i - half):i + half + 1].mean()
                for i in range(len(arr))
            ])
            prof[nm] = sm
        os.makedirs(os.path.dirname(os.path.abspath(args.profile)), exist_ok=True)
        with open(args.profile, "w", newline="") as fh:
            fh.write("\t".join(cols) + "\n")
            for i, rn in enumerate(ref_nums):
                vals = [f"{prof[nm][i]:.1f}" for nm in alignments]
                fh.write("\t".join([str(rn), ref_seq[i]] + vals) + "\n")

        # Most and least conserved windows, which is what the profile is for.
        for nm in alignments:
            sm = prof[nm]
            lo_i, hi_i = int(np.argmin(sm)), int(np.argmax(sm))
            print(f"  {nm:<20} least conserved window centred at ref {ref_nums[lo_i]} "
                  f"({sm[lo_i]:.0f}%), most at ref {ref_nums[hi_i]} ({sm[hi_i]:.0f}%)")

    # ---- is one region less conserved than another? two-proportion z-test
    if len(regions) >= 2:
        print(f"\nregion-vs-region identity, two-proportion z-test (no phylogenetic correction)")
        for nm, (pairs, seq, nums, ca, pl) in alignments.items():
            counts = {}
            for lab, lo, hi in regions:
                n = ok = 0
                for i, j in pairs:
                    if i is None or j is None or not (lo <= ref_nums[i] <= hi):
                        continue
                    n += 1
                    ok += ref_seq[i] == seq[j]
                counts[lab] = (ok, n)
            for a, b in itertools.combinations(counts, 2):
                k1, n1 = counts[a]
                k2, n2 = counts[b]
                if min(n1, n2) < 10:
                    continue
                p1, p2 = k1 / n1, k2 / n2
                p = (k1 + k2) / (n1 + n2)
                se = np.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
                if se == 0:
                    continue
                z = (p1 - p2) / se
                pv = 2 * (1 - norm_cdf(abs(z)))
                flag = "  *" if pv < 0.05 else ""
                print(f"  {nm:<18}{a:>12} {100*p1:5.1f}% (n={n1:3d})  vs {b:>12} "
                      f"{100*p2:5.1f}% (n={n2:3d})   z={z:+5.2f}  p={pv:.3f}{flag}")
    print()


def norm_cdf(x: float) -> float:
    """Standard normal CDF via erf, so scipy is not needed."""
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


if __name__ == "__main__":
    main()
