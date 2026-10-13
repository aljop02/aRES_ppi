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

Histidine is reported but never counted into the net, because its charge at
physiological pH depends on local environment rather than on composition.

Dependencies: none.

Usage
-----
  python3 charge_composition.py seq.fasta
  python3 charge_composition.py seq.fasta --regions 1-48:Zn-ribbon 49-120:NTD-1
  python3 charge_composition.py seq.fasta --window 15 --tsv charge_window.tsv
"""

from __future__ import annotations

import argparse
import sys

BASIC = "KR"
ACIDIC = "DE"

DEFAULT_REGIONS = [
    "1-48:Zn ribbon",
    "49-120:NTD-1",
    "1-120:N-block",
    "126-447:RES side",
]


def read_fasta(path: str) -> str:
    seq = []
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line)
    return "".join(seq).upper()


def parse_region(spec: str):
    """'49-120:NTD-1' -> (49, 120, 'NTD-1')"""
    rng, _, label = spec.partition(":")
    try:
        lo, hi = (int(x) for x in rng.split("-"))
    except ValueError:
        raise SystemExit(f"bad region '{spec}', expected LO-HI[:NAME]")
    return lo, hi, (label or rng)


def main() -> None:
    ap = argparse.ArgumentParser(description="Per-region charge composition.")
    ap.add_argument("fasta", help="FASTA file (single sequence)")
    ap.add_argument("--regions", nargs="*", default=DEFAULT_REGIONS,
                    metavar="LO-HI:NAME", help="regions to report")
    ap.add_argument("--window", type=int, default=0,
                    help="also emit a sliding-window net charge profile")
    ap.add_argument("--tsv", help="write the window profile here instead of stdout")
    args = ap.parse_args()

    seq = read_fasta(args.fasta)
    if not seq:
        raise SystemExit("no sequence found")
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
