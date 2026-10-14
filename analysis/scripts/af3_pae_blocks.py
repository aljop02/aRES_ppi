#!/usr/bin/env python3
"""
af3_pae_blocks.py -- block-averaged PAE from AlphaFold 3 full_data JSONs.

Why this exists: pLDDT tells you whether a region is locally well modelled, and it
will happily be high for two domains whose *relative* placement is a guess. The
predicted aligned error does tell you that, and it is the only number that licenses
a statement about inter-domain arrangement.

Concretely: superposing two orthologs on one domain and measuring another domain in
that frame can give a large RMSD either because the arrangement genuinely differs or
because neither model determined the arrangement in the first place. Block PAE
separates those two cases, so run this before believing any domain-rearrangement
claim derived from a superposition.

Rough reading of mean inter-block PAE:

  < 5 A    relative placement well determined
  5-10 A   determined in outline
  10-15 A  weak
  > 15 A   not determined -- do not interpret the relative orientation at all

Ligand tokens (a Zn ion is one token) are included, so the ion's PAE against each
protein block is reported too: that is the honest check on whether the metal's
placement is pinned by the site or merely dropped nearby.

Dependencies: numpy.

Usage
-----
  python3 af3_pae_blocks.py <dir-with-full_data-jsons> [...] \
      --blocks 1-48:Zn-ribbon 49-120:NTD1 126-447:RES 448-459:tail
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np


def parse_blocks(specs):
    out = []
    for spec in specs:
        rng, _, label = spec.partition(":")
        try:
            lo, hi = (int(x) for x in rng.split("-"))
        except ValueError:
            print(f"bad block '{spec}', expected LO-HI[:NAME]", file=sys.stderr)
            continue
        out.append((label or rng, lo, hi))
    return out


def block_masks(res_ids, chain_ids, blocks, protein_chain="A"):
    res = np.asarray(res_ids)
    ch = np.asarray(chain_ids)
    masks = {}
    for lab, lo, hi in blocks:
        masks[lab] = (ch == protein_chain) & (res >= lo) & (res <= hi)
    masks["<ligand>"] = ch != protein_chain
    return {k: v for k, v in masks.items() if v.any()}


def analyse(d: str, blocks, quiet=False):
    files = sorted(glob.glob(os.path.join(d, "*full_data_*.json")))
    if not files:
        print(f"{d}: no *full_data_*.json found", file=sys.stderr)
        return None
    name = os.path.basename(os.path.normpath(d))
    print("=" * 78)
    print(f"{name}   ({len(files)} models)")
    print("=" * 78)

    acc = defaultdict(list)
    labels = None
    for path in files:
        with open(path) as fh:
            data = json.load(fh)
        pae = np.asarray(data["pae"], dtype=float)
        masks = block_masks(data["token_res_ids"], data["token_chain_ids"], blocks)
        labels = list(masks)
        for a in labels:
            for b in labels:
                sub = pae[np.ix_(masks[a], masks[b])]
                if sub.size:
                    # Symmetrise: PAE is not symmetric, and the mean of both
                    # directions is what a "relative placement" claim rests on.
                    acc[(a, b)].append(float(sub.mean()))

    n = len(files)
    print(f"\n  mean PAE (A), averaged over {n} models, symmetrised")
    w = max(11, max(len(l) for l in labels) + 1)
    print("  " + " " * w + "".join(f"{l:>{w}}" for l in labels))
    for a in labels:
        row = f"  {a:<{w}}"
        for b in labels:
            v = np.mean(acc[(a, b)] + acc[(b, a)]) if (a, b) in acc else float("nan")
            row += f"{v:>{w}.1f}"
        print(row)

    print(f"\n  inter-block placement, worst-to-best")
    pairs = []
    for i, a in enumerate(labels):
        for b in labels[i + 1:]:
            if (a, b) not in acc:
                continue
            vals = np.array(acc[(a, b)] + acc[(b, a)])
            pairs.append((float(vals.mean()), float(vals.std()), a, b))
    for m, sd, a, b in sorted(pairs, reverse=True):
        if m > 15:
            verdict = "NOT DETERMINED -- do not interpret the relative orientation"
        elif m > 10:
            verdict = "weak"
        elif m > 5:
            verdict = "determined in outline"
        else:
            verdict = "well determined"
        print(f"    {a:>12} <-> {b:<12} {m:5.1f} A (sd {sd:.1f} over models)   {verdict}")
    print()
    return {"name": name, "pairs": pairs}


def main() -> None:
    ap = argparse.ArgumentParser(description="Block-averaged PAE from AF3 full_data JSONs.")
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--blocks", nargs="+", required=True, metavar="LO-HI:NAME")
    args = ap.parse_args()

    blocks = parse_blocks(args.blocks)
    for d in args.dirs:
        analyse(d, blocks)


if __name__ == "__main__":
    main()
