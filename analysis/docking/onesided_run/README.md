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

**D106 at 0.6% cannot be in contact with aRES.** It is not solvent-reachable at all, so it
is physically incapable of forming an intermolecular contact. Restraining it would drag the
docking somewhere wrong.

This is independent confirmation of the two-step model in `docs/05_status_after_vault_integration.md`
§2, arrived at from a completely different measurement:

> **D106 is the scaffold, S122 is the contact.** D106 caps the G121 backbone amide *inside*
> gp5 (2.71 Å in the WT model) and holds the loop in the conformation that presents S122.
> It escapes indirectly. S122, at 18.2% exposure, is partly reachable and donates its OG to
> aRES D104 in the cofold — it escapes directly.

Burial was derived here from geometry alone, with no reference to the H-bond analysis that
produced the same conclusion. Two independent routes to the same mechanism.

S122 is therefore **passive**, not active: passive residues can contribute to the interface
without being forced into it, which is the right treatment for a partly buried contact residue.

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
