#!/usr/bin/env python3
"""
homolog_compare.py -- map residue numbering and superpose across aRES orthologs.

The problem this solves: every residue number in this project (C8/C11/C31/C34,
D90, E89, E93) is in the *UC5* frame, and orthologs have different lengths --
Klebsiella aRES is 443 aa against UC5's 459. Quoting "E93" for a Klebsiella model
is meaningless until the frames are aligned, and eyeballing an alignment is how
frame errors get into documents (see docs/06 on the UC5/S_34AVR mix-up).

So: global pairwise alignment (Needleman-Wunsch with affine gaps, BLOSUM62) between
the reference and each ortholog, then

  * per-region percent identity
  * what each reference position of interest maps to in every ortholog
  * CA superposition (Kabsch) over a chosen reference region, using only
    alignment-equivalent pairs, with RMSD

Inputs are structures (.cif/.pdb) or FASTA (.fasta/.fa/.faa), mixed freely, and the
first entry is the reference frame. A structure's sequence is read from its own
coordinates, so the alignment and the coordinates cannot drift apart. A FASTA file
may hold many records, which is how a whole ortholog panel is passed in one
argument -- those entries have no coordinates, so the superposition and pLDDT
columns skip them.

Caveat on the RMSD: it is computed over aligned pairs after a least-squares fit, so
it answers "is this the same fold in the same conformation" and nothing more. A low
RMSD over a region says the models agree there; it does not make either one right.

Caveat on the identity: --identity-tsv reports it two ways. `id_aligned` divides by
the aligned columns only, so a deletion inside a region does not count against that
region; `id_allpos` divides by every reference position in the region, so it does.
They differ exactly where indels fall, and docs/10 section 2 found the N-block's
indels concentrated in NTD-1 -- so `id_aligned` is the conservative choice when the
question is whether NTD-1 diverges fastest. Both columns are written; say which one
a number came from.

Statistics: the per-ortholog two-proportion z-test at the end of this script treats
alignment positions as independent and orthologs as independent of each other.
Neither holds. For a panel, run ortholog_panel_stats.py on the --identity-tsv
output -- it tests at the ortholog level, corrects for the phylogeny, and corrects
for the several region pairs being compared.

Dependencies: numpy.

Usage
-----
  python3 homolog_compare.py ref.cif other1.cif other2.cif \
      --names UC5 E.chengduensis Klebsiella \
      --regions 1-48:Zn-ribbon 49-120:NTD1 126-447:RES \
      --map 8 11 31 34 89 90 93 \
      --super 1-120

  python3 homolog_compare.py data/sequences/ares_uc5.fasta panel.fasta \
      --regions 1-48:Zn-ribbon 49-120:NTD1 126-447:RES \
      --identity-tsv results/panel_region_identity.tsv
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


# BLOSUM62 widened by one row and column for "anything not one of the twenty",
# which scores -1 against everything including itself. Encoding the sequences
# against this lets a whole alignment row be scored in one array lookup.
_UNK = len(AA)
EXT62 = np.full((_UNK + 1, _UNK + 1), -1, dtype=np.int64)
EXT62[:_UNK, :_UNK] = BLOSUM62
_CODE = np.full(128, _UNK, dtype=np.int64)
for _i, _c in enumerate(AA):
    _CODE[ord(_c)] = _i


def encode(s: str) -> np.ndarray:
    """Sequence -> residue indices into EXT62."""
    return _CODE[np.frombuffer(s.encode("ascii", "replace"), dtype=np.uint8)]


# ------------------------------------------------------------------ alignment


def align(s1: str, s2: str):
    """Needleman-Wunsch with affine gaps (Gotoh). -> list of (i, j) 0-based pairs.

    Gaps are represented by None on one side.

    Filled a row at a time rather than a cell at a time, which matters once a
    panel is being aligned instead of two orthologs: 18 sequences against the
    reference plus the all-against-all matrix the phylogenetic correction needs is
    ~170 alignments, and the cell-at-a-time version took about a second each.

    M and X at row i depend only on row i-1, so they vectorise directly. Y depends
    on its own row -- Y[i,j] = max(M[i,j-1]+open, Y[i,j-1]+extend) -- which is a
    scan, but subtracting the per-column extension cost turns it into a cumulative
    maximum. Tie-breaking is kept identical to the cell-at-a-time recurrences
    (prefer the match state, then prefer opening a gap over extending one), so the
    traceback below is unchanged and the output is the same alignment, not merely
    an equally-scoring one.
    """
    n, m = len(s1), len(s2)
    NEG = -10**9
    e1, e2 = encode(s1), encode(s2)
    # M: match state, X: gap in s2 (consume s1), Y: gap in s1 (consume s2)
    M = np.full((n + 1, m + 1), NEG, dtype=np.int64)
    X = np.full((n + 1, m + 1), NEG, dtype=np.int64)
    Y = np.full((n + 1, m + 1), NEG, dtype=np.int64)
    pM = np.zeros((n + 1, m + 1), dtype=np.int8)
    pX = np.zeros((n + 1, m + 1), dtype=np.int8)
    pY = np.zeros((n + 1, m + 1), dtype=np.int8)

    M[0, 0] = 0
    rows = np.arange(1, n + 1)
    cols = np.arange(1, m + 1)
    X[1:, 0] = GAP_OPEN + (rows - 1) * GAP_EXTEND
    pX[1:, 0] = 1
    Y[0, 1:] = GAP_OPEN + (cols - 1) * GAP_EXTEND
    pY[0, 1:] = 2

    jv = np.arange(m + 1, dtype=np.int64)
    for i in range(1, n + 1):
        srow = EXT62[e1[i - 1], e2]

        cand = np.stack((M[i - 1, :-1], X[i - 1, :-1], Y[i - 1, :-1]))
        k = cand.argmax(0)
        M[i, 1:] = cand[k, np.arange(m)] + srow
        pM[i, 1:] = k

        c_open = M[i - 1, :] + GAP_OPEN
        c_ext = X[i - 1, :] + GAP_EXTEND
        opened = c_open >= c_ext
        X[i, :] = np.where(opened, c_open, c_ext)
        pX[i, :] = np.where(opened, 0, 1)
        X[i, 0] = GAP_OPEN + (i - 1) * GAP_EXTEND
        pX[i, 0] = 1

        # enter[j] is the score of starting the Y run at column j; removing the
        # j * extend term makes a running maximum give every Y[i, j] at once.
        enter = np.empty(m + 1, dtype=np.int64)
        enter[0] = Y[i, 0]
        enter[1:] = M[i, :-1] + GAP_OPEN
        shifted = enter - jv * GAP_EXTEND
        run = np.maximum.accumulate(shifted)
        Y[i, :] = run + jv * GAP_EXTEND
        pY[i, 1:] = np.where(shifted[1:] >= run[:-1], 0, 2)

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


FASTA_EXT = {".fasta", ".fa", ".faa", ".fas", ".seq"}


def read_fasta_records(path: str):
    """-> [(name, seq), ...]. The name is the first whitespace-delimited token of
    the header, so '>A0A2M9WDD2 Pantoea rodasii | ...' is called A0A2M9WDD2."""
    recs, name, buf = [], None, []
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = (line[1:].split() or ["unnamed"])[0], []
            elif line:
                buf.append(line.upper())
    if name is not None:
        recs.append((name, "".join(buf)))
    if not recs:
        raise SystemExit(f"{path}: no FASTA records")
    return recs


def load_inputs(paths, names):
    """-> [(name, seq, resnums, ca, plddt_by_resnum), ...], reference first.

    One path can contribute several entries, because a FASTA file may hold a whole
    ortholog panel. FASTA entries get consecutive numbering from 1 and no
    coordinates, so a NaN CA array stands in -- the superposition sections drop
    NaN pairs already, which makes them skip FASTA entries without a special case.
    """
    entries = []
    for p in paths:
        if os.path.splitext(p)[1].lower() in FASTA_EXT:
            for rec_name, seq in read_fasta_records(p):
                entries.append((rec_name, seq, list(range(1, len(seq) + 1)),
                                np.full((len(seq), 3), np.nan), {}))
        else:
            seq, nums, ca, pl = read_chain(p)
            entries.append((os.path.basename(p), seq, nums, ca, pl))
    # --names overrides in order, and may cover only the first few entries.
    out = [((names[k] if k < len(names) else e[0]),) + e[1:]
           for k, e in enumerate(entries)]
    # Everything downstream is keyed by name, so a duplicate would silently drop
    # an ortholog from the panel and from the statistics.
    seen = {}
    for nm, *_ in out:
        seen[nm] = seen.get(nm, 0) + 1
    dup = sorted(nm for nm, c in seen.items() if c > 1)
    if dup:
        raise SystemExit(f"duplicate input names: {', '.join(dup)} -- "
                         f"FASTA headers must be unique in their first token")
    return out


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
    ap.add_argument("inputs", nargs="+", metavar="INPUT",
                    help="structures (.cif/.pdb) and/or FASTA (.fasta/.fa/.faa, "
                         "may be multi-record); the first entry is the reference frame")
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
    ap.add_argument("--identity-tsv", metavar="PATH",
                    help="write per-ortholog per-region identity counts here, for "
                         "ortholog_panel_stats.py")
    ap.add_argument("--map-tsv", metavar="PATH",
                    help="write the --map table here as well; the printed version is "
                         "unreadable past a handful of orthologs")
    args = ap.parse_args()

    entries = load_inputs(args.inputs, list(args.names))
    names = [e[0] for e in entries]
    chains = [e[1:] for e in entries]
    ref_name, (ref_seq, ref_nums, ref_ca, ref_pl) = names[0], chains[0]

    print("=" * 78)
    print(f"reference frame: {ref_name}  ({len(ref_seq)} aa)")
    print("=" * 78)
    if len(entries) <= 12:
        for nm, (seq, _, _, _) in zip(names, chains):
            print(f"  {nm:<20} {len(seq)} aa")
    else:
        lens = [len(e[1]) for e in entries[1:]]
        print(f"  {len(entries) - 1} orthologs, {min(lens)}-{max(lens)} aa")

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

    # ---- per-region counts for the panel statistics
    # Counts, not just percentages, because the panel test needs to know how many
    # positions each number rests on -- a 48-residue region and a 322-residue one
    # carry very different sampling noise and the write-up has to say so.
    if args.identity_tsv:
        rows = []
        for nm, (pairs, seq, nums, ca, pl) in alignments.items():
            eq = [(i, j) for i, j in pairs if i is not None and j is not None]
            for lab, lo, hi in [("overall", None, None)] + list(regions):
                in_region = [
                    (i, j) for i, j in eq
                    if lo is None or lo <= ref_nums[i] <= hi
                ]
                n_al = len(in_region)
                n_id = sum(1 for i, j in in_region if ref_seq[i] == seq[j])
                n_ref = (len(ref_seq) if lo is None
                         else sum(1 for rn in ref_nums if lo <= rn <= hi))
                rows.append(dict(
                    ortholog=nm, region=lab,
                    lo=("" if lo is None else lo), hi=("" if hi is None else hi),
                    n_ref=n_ref, n_aligned=n_al, n_identical=n_id,
                    n_gap=n_ref - n_al,
                    id_aligned=(100.0 * n_id / n_al) if n_al else float("nan"),
                    id_allpos=(100.0 * n_id / n_ref) if n_ref else float("nan"),
                ))
        cols = ["ortholog", "region", "lo", "hi", "n_ref", "n_aligned",
                "n_identical", "n_gap", "id_aligned", "id_allpos"]
        d = os.path.dirname(os.path.abspath(args.identity_tsv))
        os.makedirs(d, exist_ok=True)
        with open(args.identity_tsv, "w", newline="") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(
                    f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c])
                    for c in cols) + "\n")
        print(f"\nper-region identity counts written to {args.identity_tsv} "
              f"({len(rows)} rows)")

    # ---- map positions of interest, machine-readable
    if args.map and args.map_tsv:
        rows = []
        for nm, (pairs, seq, nums, ca, pl) in alignments.items():
            eq = {i: j for i, j in pairs if i is not None and j is not None}
            for rn in args.map:
                try:
                    i = ref_nums.index(rn)
                except ValueError:
                    continue
                j = eq.get(i)
                rows.append(dict(
                    ortholog=nm, ref_resnum=rn, ref_aa=ref_seq[i],
                    aa=("-" if j is None else seq[j]),
                    resnum=("" if j is None else nums[j]),
                    identical=("" if j is None else int(seq[j] == ref_seq[i])),
                    plddt=("" if j is None or nums[j] not in pl
                           else f"{pl[nums[j]]:.1f}"),
                ))
        cols = ["ortholog", "ref_resnum", "ref_aa", "aa", "resnum", "identical",
                "plddt"]
        os.makedirs(os.path.dirname(os.path.abspath(args.map_tsv)), exist_ok=True)
        with open(args.map_tsv, "w", newline="") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
        print(f"\nmapped positions written to {args.map_tsv} ({len(rows)} rows)")

    # ---- map positions of interest
    # Printed unless a panel is being mapped and the TSV already has it, since the
    # printed table is one column per ortholog and stops being readable past a few.
    if args.map and (len(alignments) <= 8 or not args.map_tsv):
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
                    conf = f" pLDDT {pl[nums[j]]:.0f}" if nums[j] in pl else ""
                    cells.append(f"{seq[j]}{nums[j]}{same}{conf}")
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
    # A FASTA entry has no coordinates, so the two structural sections have nothing
    # to say about it. Both loops below skip such entries rather than printing a
    # column of "--" for every member of a panel.
    with_coords = {nm: v for nm, v in alignments.items()
                   if not np.isnan(v[3]).all()}
    if np.isnan(ref_ca).all():
        with_coords = {}
    skipped = len(alignments) - len(with_coords)
    print(f"\nCA RMSD to {ref_name} after Kabsch superposition (aligned pairs only)")
    if skipped:
        print(f"  ({skipped} sequence-only input(s) skipped: no coordinates)")
    print(f"  {'ortholog':<20}{'region':<16}{'n_pairs':>9}{'RMSD':>9}")
    for nm, (pairs, seq, nums, ca, pl) in with_coords.items():
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
    if len(targets) >= 2 and with_coords:
        print(f"\ninter-domain arrangement: fit on one region, measure another in that frame")
        print(f"  {'ortholog':<18}{'fit on':<14}{'measured':<14}{'n':>5}{'RMSD':>9}")
        for nm, (pairs, seq, nums, ca, pl) in with_coords.items():
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
