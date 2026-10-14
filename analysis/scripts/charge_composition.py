#!/usr/bin/env python3
"""
charge_composition.py -- per-region charge composition of a protein sequence.

Counts Lys+Arg against Asp+Glu over named regions and reports the net charge
both absolutely and per 100 residues, which is the number worth comparing
between regions of different length.

This was written for the aRES N-block, where the two N-terminal modules turn
out to carry almost exactly equal and opposite charge density:

    Zn ribbon 1-48    net +11   (+22.9 per 100)   ten Lys, no Glu
    NTD-1 49-120      net -16   (-22.2 per 100)   nineteen acidic, no Lys
    RES side 126-447  net  -3   ( -0.9 per 100)   an ordinary protein surface

Those numbers are UC5's, and only part of the pattern generalises. On 190 orthologs
(docs/11 section 3):

  * NTD-1's *acidity* is a family feature, and a strong one -- net-negative in
    187/190 and never net-positive in any of them, at -16.5 charges per 100
    residues against -1.2 for the PF08808 core. This is the sharpest family-level
    signal the project has; it is sharper than sequence identity, which cannot
    even separate NTD-1 from the ribbon (docs/11 section 1).
  * The equal-and-opposite *split* is a UC5 peculiarity. The ribbon's +11 is
    -5.9 on average across the family and spans -27 to +11. Do not quote the
    dipole as a property of the family.
  * "Zero Lys and exactly eleven Asp", which docs/10 section 3 reported from three
    sequences, does *not* hold: 64/190 have no Lys (a real enrichment -- 2/190 for
    the RES side -- but not an invariant), and Asp runs 3-17. Quote the charge
    *density*, not the counts.

The lesson that generalises past this protein: a count repeated in three sequences
is not a family feature, and the cost of finding that out was one run of this
script over a panel.

Histidine is reported but never counted into the net, because its charge at
physiological pH depends on local environment rather than on composition.

A FASTA with several records is treated as a panel: every record is reported
separately and then summarised across the panel, which is how the "zero Lys,
eleven Asp" claim above gets tested on more than three sequences. With --ref, the
region bounds are read in the *reference's* numbering and mapped onto each record
through a pairwise alignment, because an ortholog with a 14-residue insertion at
UC5 ref~66 does not have NTD-1 at 49-120 in its own frame, and slicing 49-120 out
of its sequence would count the wrong residues. Composition per 100 residues is
the comparable number once the mapped regions differ in length; the absolute
counts are printed next to the region length so the difference is visible.

Dependencies: numpy only when --ref is used (for the alignment), otherwise none.

Usage
-----
  python3 charge_composition.py seq.fasta
  python3 charge_composition.py seq.fasta --regions 1-48:Zn-ribbon 49-120:NTD-1
  python3 charge_composition.py seq.fasta --window 15 --tsv charge_window.tsv

  python3 charge_composition.py panel.fasta --ref data/sequences/ares_uc5.fasta \
      --regions 1-48:Zn-ribbon 49-120:NTD-1 126-447:RES \
      --tsv results/panel_charge.tsv
"""

from __future__ import annotations

import argparse
import os
import sys

BASIC = "KR"
ACIDIC = "DE"

DEFAULT_REGIONS = [
    "1-48:Zn ribbon",
    "49-120:NTD-1",
    "1-120:N-block",
    "126-447:RES side",
]


def read_fasta_records(path: str):
    """-> [(name, seq), ...]; the name is the header's first token."""
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
        raise SystemExit(f"{path}: no sequence found")
    return recs


def read_fasta(path: str) -> str:
    """Single sequence. Refuses a multi-record file rather than concatenating it,
    which would silently glue a whole panel into one 8000-residue 'protein'."""
    recs = read_fasta_records(path)
    if len(recs) > 1:
        raise SystemExit(f"{path}: {len(recs)} records -- pass --ref to treat this "
                         f"as a panel, or split the file")
    return recs[0][1]


def map_regions(ref_seq: str, seq: str, regions):
    """Translate (lo, hi, label) from the reference frame into this sequence's own
    numbering, via global pairwise alignment. -> [(lo, hi, label, n_ref_aligned)].

    The mapped bounds are the first and last aligned positions of the region, so an
    insertion inside the region widens it -- which is the point: those inserted
    residues are part of that module in this ortholog and their charge counts.
    """
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from homolog_compare import align  # numpy; only needed on this path

    pairs = align(ref_seq, seq)
    eq = {i: j for i, j in pairs if i is not None and j is not None}
    out = []
    for lo, hi, label in regions:
        js = [eq[i] for i in range(lo - 1, min(hi, len(ref_seq))) if i in eq]
        if not js:
            out.append((None, None, label, 0))
            continue
        out.append((min(js) + 1, max(js) + 1, label, len(js)))
    return out


def counts(s: str):
    k, r, h = s.count("K"), s.count("R"), s.count("H")
    d, e = s.count("D"), s.count("E")
    return dict(len=len(s), K=k, R=r, H=h, D=d, E=e,
                basic=k + r, acidic=d + e, net=k + r - d - e)


def parse_region(spec: str):
    """'49-120:NTD-1' -> (49, 120, 'NTD-1')"""
    rng, _, label = spec.partition(":")
    try:
        lo, hi = (int(x) for x in rng.split("-"))
    except ValueError:
        raise SystemExit(f"bad region '{spec}', expected LO-HI[:NAME]")
    return lo, hi, (label or rng)


def panel_report(args, records, regions) -> None:
    """Per-region composition for every record, then the across-panel summary.

    The summary is the part that tests a family-level claim: a count that is
    identical in three sequences at 34-69% identity (docs/10 section 3) is either
    a property under selection or an artifact of a small sample, and the only way
    to tell is to widen the sample and report the whole distribution rather than
    the modal value.
    """
    ref_seq = None
    if args.ref:
        ref_recs = read_fasta_records(args.ref)
        if len(ref_recs) != 1:
            raise SystemExit(f"{args.ref}: reference must hold exactly one sequence")
        ref_name, ref_seq = ref_recs[0]
        print(f"\nreference frame: {ref_name} ({len(ref_seq)} aa); region bounds "
              f"mapped onto each record by alignment")
    print(f"{args.fasta}: {len(records)} record(s)\n")

    rows = []
    for name, seq in records:
        mapped = (map_regions(ref_seq, seq, regions) if ref_seq is not None
                  else [(lo, hi, lab, hi - lo + 1) for lo, hi, lab in regions])
        for lo, hi, label, n_al in mapped:
            if lo is None:
                rows.append(dict(record=name, region=label, lo="", hi="",
                                 n_ref_aligned=0, len=0))
                continue
            c = counts(seq[lo - 1:hi])
            rows.append(dict(record=name, region=label, lo=lo, hi=hi,
                             n_ref_aligned=n_al, **c))

    labels = [lab for _, _, lab in regions]
    hdr = (f"  {'record':<14}{'region':<14}{'bounds':>12}{'len':>5}{'K':>4}{'R':>4}"
           f"{'H':>4}{'D':>5}{'E':>4}{'net':>6}{'/100':>8}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for r in rows:
        if not r["len"]:
            print(f"  {r['record']:<14}{r['region']:<14}{'-- unaligned --':>12}")
            continue
        bounds = "{}-{}".format(r["lo"], r["hi"])
        print(f"  {r['record']:<14}{r['region']:<14}"
              f"{bounds:>12}{r['len']:>5}{r['K']:>4}{r['R']:>4}"
              f"{r['H']:>4}{r['D']:>5}{r['E']:>4}{r['net']:>+6}"
              f"{r['net'] / r['len'] * 100:>+8.1f}")

    print(f"\nacross the panel, by region ({len(records)} records)")
    print(f"  {'region':<14}{'len':>9}{'K':>12}{'R':>12}{'D':>12}{'E':>12}"
          f"{'net/100':>16}")
    for lab in labels:
        sel = [r for r in rows if r["region"] == lab and r["len"]]
        if not sel:
            continue

        def rng(key):
            v = [r[key] for r in sel]
            return f"{min(v)}-{max(v)}" if min(v) != max(v) else f"{min(v)}"

        per100 = [r["net"] / r["len"] * 100 for r in sel]
        print(f"  {lab:<14}{rng('len'):>9}{rng('K'):>12}{rng('R'):>12}"
              f"{rng('D'):>12}{rng('E'):>12}"
              f"{f'{min(per100):+.1f} to {max(per100):+.1f}':>16}")

    print(f"\n  how many records have zero of a residue, by region")
    print(f"  {'region':<14}{'zero K':>10}{'zero R':>10}{'zero D':>10}{'zero E':>10}"
          f"{'net < 0':>10}")
    for lab in labels:
        sel = [r for r in rows if r["region"] == lab and r["len"]]
        if not sel:
            continue
        n = len(sel)
        cells = [sum(1 for r in sel if r[x] == 0) for x in ("K", "R", "D", "E")]
        neg = sum(1 for r in sel if r["net"] < 0)
        print(f"  {lab:<14}" + "".join(f"{f'{c}/{n}':>10}" for c in cells)
              + f"{f'{neg}/{n}':>10}")
    print()

    if args.tsv:
        cols = ["record", "region", "lo", "hi", "n_ref_aligned", "len",
                "K", "R", "H", "D", "E", "basic", "acidic", "net", "net_per100"]
        d = os.path.dirname(os.path.abspath(args.tsv))
        os.makedirs(d, exist_ok=True)
        with open(args.tsv, "w", newline="") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                r = dict(r)
                r.setdefault("K", "")
                r["net_per100"] = (f"{r['net'] / r['len'] * 100:.2f}"
                                   if r.get("len") else "")
                fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
        print(f"per-record counts written to {args.tsv}\n")


def main() -> None:
    ap = argparse.ArgumentParser(description="Per-region charge composition.")
    ap.add_argument("fasta", help="FASTA file; several records are treated as a panel")
    ap.add_argument("--regions", nargs="*", default=DEFAULT_REGIONS,
                    metavar="LO-HI:NAME", help="regions to report")
    ap.add_argument("--ref", metavar="FASTA",
                    help="reference sequence whose numbering the regions are in; "
                         "each record's bounds are mapped onto it by alignment")
    ap.add_argument("--window", type=int, default=0,
                    help="also emit a sliding-window net charge profile")
    ap.add_argument("--tsv", help="write per-record per-region counts here (or the "
                                  "window profile, with --window)")
    args = ap.parse_args()

    regions = [parse_region(s) for s in args.regions]
    records = read_fasta_records(args.fasta)

    if args.ref or len(records) > 1:
        panel_report(args, records, regions)
        return

    seq = records[0][1]
    print(f"\n{args.fasta}: {len(seq)} residues\n")

    hdr = (f"{'region':<20}{'len':>5}{'K':>4}{'R':>4}{'H':>4}"
           f"{'D':>5}{'E':>4}{'basic':>7}{'acidic':>7}{'net':>6}{'/100':>8}")
    print(hdr)
    print("-" * len(hdr))
    for spec in args.regions:
        lo, hi, label = parse_region(spec)
        s = seq[lo - 1:hi]
        if not s:
            print(f"{label:<20}  -- no residues in range")
            continue
        k, r, h = s.count("K"), s.count("R"), s.count("H")
        d, e = s.count("D"), s.count("E")
        net = k + r - d - e
        print(f"{label:<20}{len(s):>5}{k:>4}{r:>4}{h:>4}{d:>5}{e:>4}"
              f"{k + r:>7}{d + e:>7}{net:>+6}{net / len(s) * 100:>+8.1f}")
    print()

    if args.window:
        w = args.window
        rows = ["pos\tnet"]
        half = w // 2
        for i in range(len(seq)):
            s = seq[max(0, i - half):i + half + 1]
            rows.append(f"{i + 1}\t{sum(c in BASIC for c in s) - sum(c in ACIDIC for c in s)}")
        out = "\n".join(rows)
        if args.tsv:
            with open(args.tsv, "w") as fh:
                fh.write(out + "\n")
            print(f"window profile (w={w}) written to {args.tsv}\n")
        else:
            print(out)


if __name__ == "__main__":
    main()
