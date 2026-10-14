# Blind ab-initio run **with Zn modelled** — HADDOCK job 780832

The blind run from [`../blind_run/`](../blind_run/README.md) repeated with the Zn(Cys)₄ site
occupied. This is the run that answers point 3 of
[`docs/07`](../../../docs/07_docking_run_review.md) §1 — "there is no zinc in the docking run …
if docking contacts land on residues 1–48, discount them until the run is repeated with Zn
modelled."

Results and interpretation: **[`docs/09_blind_zn_run_results.md`](../../../docs/09_blind_zn_run_results.md)**.

## What is here

| | |
|---|---|
| `cluster<N>_<1-4>.pdb` | the 28 cluster representatives (top 4 of each of 7 clusters), from `780832-new_aRES_summary.tgz` — **gitignored**, per the repo policy of keeping summaries rather than raw runs |
| `plots/*.html` | the server's Bokeh analysis plots — **gitignored**; the numbers inside them are extracted to `results/blind_zn_cluster_stats.tsv` |

The downloaded archive is a **plain tar despite the `.tgz` name** (the browser gunzipped it in
transit), so `tar -xzf` fails on it and `tar -xf` is what works.

## System as docked

| | segid | contents |
|---|---|---|
| aRES (UC5) | **A** | residues **1–447** (disordered 448–459 trimmed, as recommended) **+ Zn²⁺ as residue 500, `ZN2`** |
| T7 gp5 | **B** | residues 1–704 |

The zinc is an `ATOM` record, not `HETATM`, with atom and residue name `ZN+2`/`ZN2` — so a
naive `grep HETATM` finds nothing and it is easy to conclude the metal is absent. It is present.

### The metal site is correctly built

Verified on `cluster1_1.pdb` (identical in all 28):

| | |
|---|---|
| Zn²⁺ coordinated by | **Cys8 2.17 Å, Cys11 2.08 Å, Cys31 2.27 Å, Cys34 2.12 Å** |
| next nearest sulfur | Cys78 at **14.3 Å** — nothing else is close |
| mean SG···SG edge | **3.53 Å** (was 3.67 Å in the apo AlphaFold model; ideal 3.76 Å) |

Exactly the four cysteines predicted by `zn_site_check.py`, no spurious fifth ligand, and the
tetrahedron contracted onto the metal as it should. The site is right.

**The Zn²⁺ makes no contact with gp5 in any of the 28 models** (0.0% of interface weight). Its
effect on the docking is through the aRES fold and the ribbon's electrostatics, not through any
direct metal–gp5 interaction — which also means no cluster is being held together by the ion.

## Reproducing the analysis

```bash
# per-cluster domain breakdown, with the N-terminal fray separated from the ribbon fold
python3 ../../scripts/interface_contacts.py cluster6_*.pdb --chains A B \
    --domain A:1-5:Nterm_fray --domain A:6-48:Zn_ribbon --domain A:49-120:NTD1 \
    --domain A:126-447:RES --roi B:100-125:epitope_loop

# cluster populations and whether the best structure is clustered (NOT in the summary tgz)
python3 ../../scripts/haddock_cluster_stats.py plots/fcc_hs.html \
    --tsv ../../../results/blind_zn_cluster_stats.tsv
```

`--domain A:1-5:Nterm_fray` matters. Residues 1–5 carry the charged free N-terminus, and in two
clusters they supply ~27% of the aRES interface — splitting them out is what separates a real
N-block pose from a terminus artifact.

## ⚠️ The best-scoring structure is not in this archive

The summary download contains cluster representatives only. The lowest-energy structure of the
whole run has **HADDOCK score −154.0** — better than any clustered model — and it is
**unclustered**, so it is absent here. Fetch it from the server (`structures/it1/water/`) before
the two-week deletion; it is the single most informative missing file.
