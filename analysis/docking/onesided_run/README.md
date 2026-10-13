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

## ⚠️ Status: NOT YET SUBMITTED (as of 13-Oct-2026)

The two runs that came back on 13-Oct-2026 are both **blind** runs (`aRES-T7` = job 779823 and
`1_aRES_T7` = job 780456); see
[`docs/09_haddock_blind_run_results.md`](../../../docs/09_haddock_blind_run_results.md) §6. Neither
is `2_ares_T7`. Both carry `Eair = 0` in every model, which `ranair: false` with real AIRs would not.

So this run is still outstanding, and it is now the one that matters: the blind runs failed to
converge, which is exactly the outcome `docs/07` §3 predicted and the reason this design exists.

**When it returns, check `air_boxplot.html` first.** If it says *"No Eair was calculated for this
run"*, the restraints did not take effect and it is another blind run whatever was submitted.
