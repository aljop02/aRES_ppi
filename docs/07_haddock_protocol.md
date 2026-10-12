# What to run tonight

Purely computational work, in priority order. The HADDOCK run is §2; §1 is cheaper, runs in
parallel, and is arguably more important.

---

## 1. First, in parallel and free: AF3 seed replicates

**Everything structural in this project currently rests on one cofold.** Before investing a night of
docking in that interface, find out whether it is reproducible.

Run **5 AF3 jobs of aRES(UC5) + T7 gp5 with different seeds** on the
[AlphaFold Server](https://alphafoldserver.com/), and the same for **aRES + gp5 + TrxA**. For each:

- mean ± SD ipTM and pTM
- is the **S122(gp5) ↔ D104(aRES)** contact present in all 5, some, or one?
- interface PAE around gp5 100–125
- pLDDT of gp5 100–125 and of the aRES N-block

This is MD-plan step 1 and it is the right instinct. The answer changes how much weight everything
else can carry:

| Outcome | What it means |
|---|---|
| Contact in 5/5, good interface PAE | The cofold is a real prediction. Dock and simulate with confidence. |
| Contact in 2–3/5 | Plausible but soft. The one-sided HADDOCK run below becomes the deciding evidence. |
| Contact in 1/5 | You have one lucky model, not a prediction. Say so, and lean on the wet-lab tests. |

Record it in `analysis/afold/`. Once the models are in the repo,
`interface_contacts.py` gives the per-seed contact table in one command:

```bash
python3 analysis/scripts/interface_contacts.py analysis/afold/seeds/*.cif \
    --chains A B --domain B:1-120:NTD1 --domain B:126-447:RES \
    --roi A:100-125:loop --min-plddt 70
```

---

## 2. The HADDOCK run to start tonight

### Not the blind full × full run

The blind run already going stays as the unbiased control, but do not start a second one. aRES (459)
× gp5 (704) ab initio has a severe sampling problem — size is not the limit (the server allows
~50–75k atoms, you need ~18k), the search space is.

### Run it one-sided instead

> **Restrain only the gp5 loop. Leave aRES completely free.**

This is the best design available, because the information you inject is **experimentally justified
and independent of any prediction**: the escape mutations S122P and D106G define the gp5 epitope by
genetics. The aRES side — the actual question — stays unconstrained.

So whichever aRES domain HADDOCK puts at the interface is **a prediction, not an assumption**. And
because HADDOCK's physics-based scoring shares essentially no assumptions with AlphaFold's learned
potential, agreement between them is real evidence rather than two views of the same bias.

### Pre-register the readout — before you look

Write this down now:

| Prediction | Supported if | Falsified if |
|---|---|---|
| **NTD-1 is the binding domain** | ≥60% of aRES-side interface weight falls in **1–120** | majority falls in **126–447** (RES side) |
| **D104 is in the interface** | aRES **D104** appears in the top contacts across clusters | absent from all clusters |

A one-sided run that independently lands on NTD-1 and recovers D104 is a genuinely strong result.
One that lands on the RES domain is also a result — it would mean the cofold's NTD-1 interface is
not robust, and `docs/02` §8.1 (recognition lives in RES) becomes live.

### ⚠️ Use unbound monomers, not the cofold chains

**Do not split the cofold into two files and dock those.** They are already in the bound pose, so
side chains and the loop are pre-organized for each other. HADDOCK randomizes starting *orientation*
but not conformation, so you would be doing bound-bound docking — far easier than the real problem,
and the result would not mean what you want it to mean.

Instead predict the monomers separately:
- **aRES alone** (AF3, UC5, 1–459 or tail-trimmed 1–447)
- **gp5 alone**, or **gp5 + TrxA** as a two-chain ensemble

Then it is a true unbound docking test.

### Setup

```bash
# 1. Restraints: gp5 loop active, aRES free.
python3 analysis/scripts/make_haddock_restraints.py \
    --mol-a gp5_unbound.pdb  --seg-a A --active-a 100-125 \
    --mol-b ares_unbound.pdb --seg-b B --passive-b surface \
    --out analysis/docking/onesided/air_onesided.tbl

# 2. Sanity-check which loop residues are actually exposed before trusting the run.
python3 analysis/scripts/make_haddock_restraints.py \
    --mol-a gp5_unbound.pdb --seg-a A --report 100-125
```

The script computes Shrake-Rupley relative SASA on the **isolated** chain and drops buried residues
from the active list — restraining a buried residue pulls the docking somewhere wrong. Check the
report output: if most of 100–125 comes back buried, the AF model has the loop packed against the
body rather than exposed, and that is itself worth knowing before the run.

HADDOCK 2.4 server settings (guru interface):

| Setting | Value | Why |
|---|---|---|
| Molecule A | gp5 (or gp5+TrxA) | |
| Molecule B | aRES | |
| AIR file | `air_onesided.tbl` | from the script |
| **Random removal of restraints** | **off** | with few actives, the default 50% removal throws away most of your signal |
| Semi-flexible segments | gp5 **100–125** | the epitope loop; let it adapt |
| Fully flexible segments | gp5 100–125 | the loop is disordered in the crystal (1T7P deletes 118–123), so treat it as genuinely mobile |
| Structures for rigid-body (it0) | **10000** (default 1000) | harder-than-average case; this is what makes it an overnight run |
| it1 / water | 400 / 400 | |
| Clustering | FCC, cutoff 0.6, min 4 | |

If a local **haddock3** install is available, prefer it — you control sampling and can run the
variants below in one batch. The 2.4 server is fine otherwise; queue time is the unknown, so submit
early.

### Second run, if there is capacity

**NTD-1 (aRES 1–120) × gp5, fully blind — no restraints at all.** A 120-residue partner is a
tractable ab-initio problem in a way the 459-residue one is not. If a completely unrestrained run on
the isolated domain finds the gp5 loop, that is the strongest computational evidence obtainable here,
because nothing about the answer went in.

### Analysis the next morning

```bash
python3 analysis/scripts/interface_contacts.py analysis/docking/onesided/*.pdb \
    --chains A B \
    --domain B:1-120:NTD1 --domain B:126-447:RES \
    --roi A:100-125:loop \
    --tsv results/onesided_contacts.tsv
```

Read **all** clusters, not just cluster 1. The domain breakdown on chain B answers the pre-registered
question directly; the top-pairs list tells you whether D104 came back.

---

## 3. Also worth doing while you wait — and genuinely cheap

**T2 vs T4 gp43 alignment.** UC5 aRES discriminates ~3 logs between two closely related T-even
polymerases (`docs/06` §3.1); GN04826 discriminates 1 vs 1,000. Align the two gp43 sequences, find
the differences, and check whether they cluster in the region equivalent to T7's 100–125 loop. This
is a natural escape-mutant pair that already exists in nature and costs nothing but an alignment.
**Upload the two sequences and I can run it here.**

**Check whether aRES D104 is conserved.** If the equivalent aspartate is present in T7-responsive
aRES proteins (UC5, S_34AVR) and absent or substituted in GN04826, which does **not** respond to T7,
that is an independent evolutionary argument for the interface — obtained without any experiment.
A three-sequence alignment is enough to look; the full homolog set (Tables S2/S3) makes it rigorous.
**This is the highest-value thing on the list after §1, and it needs only three sequences.**

**Map the predicted interface onto the aRES AF model** and check that the candidate surface is on
the opposite face from the RES active site (R182/R184/E236). If the predicted gp5 interface overlaps
the NAD⁺ pocket, the autoinhibition model in `docs/02` §7 gains support and the mechanism changes.

---

## 4. Honest expectation

The one-sided run can fail to converge, and a smeared-out set of clusters is a plausible outcome.
That is not wasted time — it tells you the physics-based method cannot discriminate here, which is
itself information about how much weight the AF3 cofold can carry alone. Write it down either way.

What would make the night's work count regardless of how HADDOCK goes: **the AF3 seed replicates in
§1 and the D104 conservation check in §3.** Both are cheap, neither depends on docking converging,
and either could change the project's direction on its own.
