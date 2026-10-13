# One-sided restrained run (`2_ares_T7`)

`job_params_restrained.json`. Submit via **Submit file → Full parameter file**.

## Design

Only the gp5 epitope is restrained. **The aRES side is left completely free**, so whichever
aRES domain HADDOCK places at the interface is a prediction, not an assumption. The
information injected (the gp5 loop is the epitope) comes from escape-mutant genetics, not
from any structure prediction, so a result agreeing with the AF3 cofold is independent
support rather than a second view of the same bias.

| | segid | active | passive |
|---|---|---|---|
| T7 gp5 | B | **9** residues | 25 residues |
| aRES | A | **0 — free** | 304 (whole surface) |

```
gp5 active : 103, 104, 105, 111, 112, 113, 114, 117, 118
gp5 passive: 96, 99, 100, 101, 107, 108, 109, 110, 116, 119, 120, 121, 122, 124,
             126, 130-138, 140
```

Settings: `ranair: false`, funnel 5000 / 400 / 400 / 400, seed 2002, FCC 0.6 / min 4.
Sampling is lower than the blind run because the gp5 side is pinned, which shrinks the
search space substantially.

## ⚠️ Why neither escape residue is in the active list

Relative solvent accessibility of the loop, computed on the isolated gp5 chain:

| residue | rel. SASA | |
|---|---|---|
| **D106** | **0.6%** | escape mutation — **completely buried** |
| **S122** | **18.2%** | escape mutation — partly exposed, below the 40% active cutoff |
| K103, D104, T105, R111, S112, G113, K114, G117, K118 | 42–80% | → active |

**D106 at 0.6% cannot be in contact with aRES.** It is not solvent-reachable at all, so
restraining it would drag the docking somewhere wrong.

This is consistent with — but does **not** establish — the intramolecular role of D106. That
was already hypothesis H1 in `notes/vault/MD Plan ...md` ("Asp106 carboxylate caps the Gly121
N-H … D106G removes the cap and adds backbone flexibility"), proposed before this calculation
was run. The SASA figures were computed to decide the restraint lists, and that is what they
settle. They rule out D106 being *additionally* a direct contact, which nobody had proposed.

The more informative number is **S122 at 18.2%**. Partly exposed means the OG→aRES D104
hydrogen bond is geometrically possible but not generous — a real constraint on H2, since the
contact residue does not sit freely on the surface.

S122 is therefore **passive**, not active: passive residues can contribute to the interface
without being forced into it, which is the right treatment for a partly buried contact residue.

**What would actually test H1** is a side-chain-independent comparison across the WT, D106G and
S122P models — Cα(106)–Cα(121) distance plus backbone RMSD over 100–125 — which is defined for
glycine too, unlike the OD2/CB/N measurements in the three PyMOL panels
(`docs/05` §3). If D106G acts by destabilising the loop while S122P deletes a contact, the two
should also behave *differently* in a pull-down; if both simply fail to bind, the two-step
picture is too simple.

## Analysis when it finishes

```bash
python3 ../../scripts/interface_contacts.py <run>/structures/it1/water/*.pdb \
    --chains B A --domain A:1-120:NTD1 --domain A:126-447:RES \
    --roi B:100-125:loop --tsv ../../../results/onesided_contacts.tsv
```

**Pre-registered readout** (`docs/07_docking_run_review.md` §3) — commit to it before looking:

| | supported | falsified |
|---|---|---|
| NTD-1 is the binding domain | ≥60% of aRES interface weight in **1–120** | majority in **126–447** |
| D104 is in the interface | aRES **D104** in the top contacts across clusters | absent from all clusters |

Unlike the blind run, this one localises the aRES side whatever the outcome.

---

## ⚠️ Outcome: job 780457 FAILED — the partners never touched

Submitted and returned 13-Oct-2026 as job **780457** (`2_ares_T7`). It did not produce a docking
result at all. Full analysis in
[`docs/09_haddock_blind_run_results.md`](../../../docs/09_haddock_blind_run_results.md) §10.

**In all 400 structures the two chains sit 266–283 Å apart**, on the x-axis at y,z ≈ 0 — their
untouched starting positions. Every intermolecular term is exactly zero, buried surface area is the
`-999999` error sentinel, and the 400 structures form 400 clusters of size 1.

**Root cause: the generated AIR table was empty.** With `ranair: false`, `cmrest: false` and
`surfrest: false`, nothing else in HADDOCK's energy function pulls separated molecules together —
at 217 Å, vdW and electrostatics are far past cutoff, so there is not even a gradient to follow.

The restraint *definition* in `job_params_restrained.json` is valid: every restrained residue
exists, segids and chain IDs are consistent, `inter_mat` is all 1.0. The failure is that the server
never turned those lists into restraints. Prime suspect: **partner 1 has `activereslist = []`**.
HADDOCK builds AIRs outward from active residues, so with none on aRES only gp5's 9 should have
produced restraints — and none appeared. The "Submit file → Full parameter file" route may also skip
the restraint builder entirely.

### The fix — do not let the server infer the restraints

[`ambig.tbl`](ambig.tbl) is now generated explicitly and committed:

```bash
python3 ../../scripts/make_ambig_tbl.py \
    --from-params job_params_restrained.json --source-segid B -o ambig.tbl
```

9 restraints, one per gp5 active residue, each targeting all 304 aRES surface residues. Upload it as
the **ambiguous restraints** file when resubmitting.

### Before spending queue time again

1. **Confirm `ambig.tbl` is non-empty** and that the server accepted it (9 `assign` blocks).
2. **Check the segids.** The uploaded PDBs carry blank segid columns (73–76) and rely on the server
   filling them from the chain ID. If that ever fails, `segid A` / `segid B` match nothing and you
   get this exact failure again. Setting segids explicitly in the input PDBs removes the risk.
3. **Set `noecv: false`.** With only 9 AIRs, the default 50% random removal for cross-validation
   (`ncvpart: 2.0`) leaves very few restraints per trial.
4. **Run a 100-structure test job first.** Ten minutes, and `Eair ≠ 0` plus `BSA > 0` in the output
   proves the restraints are live before you commit 5000 structures to it.

**The check that caught this:** `Eair = 0`. Keep using it. It is the single diagnostic that
separates "docked badly" from "never docked".
