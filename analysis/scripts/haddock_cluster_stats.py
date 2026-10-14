#!/usr/bin/env python3
"""
haddock_cluster_stats.py -- recover the full per-cluster statistics of a HADDOCK 2.4
run from the Bokeh analysis plots the server ships.

Why this exists
---------------
The `*_summary.tgz` download gives you only the top four structures of each cluster.
That hides the two numbers that decide whether a run converged:

  * the cluster *populations* (a cluster of 8 out of 207 is not a result)
  * whether the best-scoring structure is clustered at all

Both are present in the Bokeh plot HTML (`fcc_hs.html` etc.), which embeds the whole
water-refined analysis set as a base64 float64 ColumnDataSource. This script decodes
it and prints/writes the per-cluster table.

HADDOCK deletes results from the server after two weeks, so run this while the plots
are still downloadable -- afterwards the populations are unrecoverable.

Note on `fcc` and `irmsd`: both are measured against the run's *reference* structure
(the lowest-energy model), not within-cluster. FCC near zero for every cluster
therefore means no cluster shares contacts with the best-scoring pose.

Dependencies: numpy only.

Usage
-----
  python3 haddock_cluster_stats.py plots/fcc_hs.html
  python3 haddock_cluster_stats.py plots/fcc_hs.html --tsv results/run_cluster_stats.tsv
"""

from __future__ import annotations

import argparse
import base64
import json
import re
import sys

import numpy as np

# HADDOCK 2.4 water-refined score.
HS_WEIGHTS = {"vdw": 1.0, "elec": 0.2, "desolv": 1.0, "air": 0.1}

NUMERIC = ("hs", "vdw", "elec", "desolv", "air", "fcc", "irmsd", "lrmsd", "rmsd")


def _decode(col):
    """Bokeh ships numeric columns as base64 float64; categorical ones as plain lists."""
    if isinstance(col, list):
        return np.array(col, dtype=object)
    if isinstance(col, dict) and "__ndarray__" in col:
        raw = base64.b64decode(col["__ndarray__"])
        return np.frombuffer(raw, dtype=col.get("dtype", "float64"))
    raise ValueError(f"unrecognised column encoding: {type(col)}")


def load_source(path: str) -> dict:
    """Pull the ColumnDataSource carrying per-structure statistics out of the HTML."""
    html = open(path, encoding="utf-8", errors="replace").read()
    blobs = re.findall(
        r'<script type="application/json"[^>]*>(.*?)</script>', html, re.S
    )
    if not blobs:
        sys.exit(f"{path}: no embedded Bokeh JSON found (is this a HADDOCK plot?)")

    doc = json.loads(blobs[0])
    hits: list[dict] = []

    def walk(obj):
        if isinstance(obj, dict):
            data = obj.get("data") or obj.get("attributes", {}).get("data")
            if isinstance(data, dict) and "clus" in data and "hs" in data:
                hits.append(data)
            for value in obj.values():
                walk(value)
        elif isinstance(obj, list):
            for value in obj:
                walk(value)

    walk(doc)
    if not hits:
        sys.exit(f"{path}: embedded JSON has no cluster/score data source")
    return hits[0]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plot", help="a HADDOCK Bokeh plot, e.g. plots/fcc_hs.html")
    ap.add_argument("--tsv", help="write the per-cluster table here")
    args = ap.parse_args()

    src = load_source(args.plot)
    cols = {k: _decode(v) for k, v in src.items()}

    clus = np.array([int(c) for c in cols["clus"]])
    legend = np.array([str(x) for x in cols["legend"]])
    n = len(clus)

    # "Other" is HADDOCK's label for structures that fell outside every cluster.
    unclustered = int((legend == "Other").sum())

    print(f"structures in the water-refined analysis set : {n}")
    print(f"clustered                                    : {n - unclustered}"
          f" ({100 * (n - unclustered) / n:.0f}%)")
    print(f"unclustered ('Other')                        : {unclustered}"
          f" ({100 * unclustered / n:.0f}%)")

    best = int(np.argmin(cols["hs"].astype(float)))
    print(f"\nbest-scoring structure : HADDOCK score {float(cols['hs'][best]):.1f}"
          f"  -> {legend[best]}")
    if legend[best] == "Other":
        print("  *** the best structure is UNCLUSTERED: the run did not converge on it.")

    rows = []
    for c in sorted(set(clus)):
        m = clus == c
        label = sorted(set(legend[m]))
        name = label[0] if len(label) == 1 else ",".join(label)
        hs = cols["hs"][m].astype(float)
        row = {
            "cluster": name,
            "n": int(m.sum()),
            "pct_of_set": 100.0 * int(m.sum()) / n,
            "hs_mean": hs.mean(),
            "hs_sd": hs.std(ddof=1) if m.sum() > 1 else 0.0,
            "hs_best": hs.min(),
        }
        for k in ("irmsd", "fcc", "elec", "vdw", "desolv", "air"):
            if k in cols:
                row[k] = cols[k][m].astype(float).mean()
        rows.append(row)
    rows.sort(key=lambda r: r["hs_mean"])

    hdr = (f"\n{'cluster':>12} {'n':>4} {'%set':>5} {'HS mean':>9} {'sd':>6} "
           f"{'HS best':>8} {'i-RMSD':>7} {'FCC':>5} {'elec':>8} {'vdw':>7} "
           f"{'desolv':>7} {'AIR':>6}")
    print(hdr)
    print("-" * (len(hdr) - 1))
    for r in rows:
        print(f"{r['cluster']:>12} {r['n']:>4} {r['pct_of_set']:>4.0f}% "
              f"{r['hs_mean']:>9.1f} {r['hs_sd']:>6.1f} {r['hs_best']:>8.1f} "
              f"{r.get('irmsd', float('nan')):>7.1f} {r.get('fcc', float('nan')):>5.2f} "
              f"{r.get('elec', float('nan')):>8.1f} {r.get('vdw', float('nan')):>7.1f} "
              f"{r.get('desolv', float('nan')):>7.2f} {r.get('air', float('nan')):>6.1f}")

    print("\ni-RMSD and FCC are measured against the run's lowest-energy structure,")
    print("not within-cluster. FCC ~ 0 across the board means the clusters share no")
    print("contacts with the best-scoring pose.")

    if args.tsv:
        keys = ["cluster", "n", "pct_of_set", "hs_mean", "hs_sd", "hs_best",
                "irmsd", "fcc", "elec", "vdw", "desolv", "air"]
        with open(args.tsv, "w", encoding="utf-8") as fh:
            fh.write("\t".join(keys) + "\n")
            for r in rows:
                fh.write("\t".join(
                    f"{r[k]:.4f}" if isinstance(r.get(k), float) else str(r.get(k, ""))
                    for k in keys) + "\n")
        print(f"\nwrote {args.tsv}")


if __name__ == "__main__":
    main()
