# Results of the blind run with Zn modelled (HADDOCK 780832)

Run record and the zinc validation: [`analysis/docking/blind_zn_run/`](../analysis/docking/blind_zn_run/README.md).

**Verdict against the pre-registered readout in [`docs/07`](07_docking_run_review.md) §3: the run
did not converge, so by the rule written down in advance the correct conclusion is _no
conclusion_.** Two features are nonetheless stable across every cluster, and one cluster is clean
enough to generate a specific, cheap experiment. Both are flagged below as hypothesis-generating,
not as support.

---

## 1. The run did not converge

The numbers that decide this are **not in the summary download** — they had to be recovered from
the server's plots (`analysis/scripts/haddock_cluster_stats.py`, output in
`results/blind_zn_cluster_stats.tsv`):

| cluster | n | % of set | HADDOCK score (mean ± sd) | best | i-RMSD | FCC |
|---|---|---|---|---|---|---|
| 5 | 9 | 4% | −101.3 ± 17.4 | −120.7 | 21.4 | 0.00 |
| 6 | 8 | 4% | −99.9 ± 9.1 | −109.4 | 17.7 | 0.00 |
| 2 | 13 | 6% | −88.0 ± 11.0 | −107.0 | 38.2 | 0.00 |
| **unclustered** | **75** | **36%** | −69.3 ± 25.7 | **−154.0** | 29.2 | 0.02 |
| 1 | 73 | 35% | −64.9 ± 15.3 | −100.0 | 23.3 | 0.00 |
| 4 | 10 | 5% | −63.6 ± 14.4 | −92.3 | 23.3 | 0.00 |
| 7 | 7 | 3% | −56.4 ± 7.4 | −69.1 | 43.2 | 0.00 |
| 3 | 12 | 6% | −47.8 ± 6.7 | −56.8 | 42.0 | 0.00 |

Three independent signs of non-convergence:

1. **The best-scoring structure in the run is unclustered**, at −154.0 — better than the best
   member of any cluster. The top pose is a singleton.
2. **i-RMSD 17.7–43.2 Å** between clusters and the reference, and **FCC ≈ 0.00 for every
   cluster**: no cluster shares contacts with the best-scoring pose. These are seven unrelated
   binding modes, not seven views of one.
3. **36% of the analysis set fell outside every cluster**, and the two best-scoring clusters hold
   9 and 8 structures out of 207.

Population and score also disagree: cluster 1 is by far the largest (73, 35%) but ranks fifth by
score, while clusters 5 and 6 score best on 4% of the set each. In a converged run these agree.

**This is the expected outcome at this size and it is not evidence against the hypothesis** — it
is an unsuccessful search, exactly as `docs/07` §3 said in advance. What follows is read in that
light.

### The score is almost entirely electrostatic

`HADDOCK score = vdw + 0.2·elec + desolv + 0.1·AIR`. For cluster 5, `0.2 × (−425.9) = −85` of a
−101 score. Cluster 5 has the best score and the *worst* van der Waals term of the N-block
clusters (−21.6 vs −45.4 for cluster 6), i.e. it is winning on long-range charge, not on packing.
With two large, patchy-charged surfaces this is the classic way a blind run converges on "the most
charged patch" rather than the biological interface. Cluster 6's score is the more trustworthy of
the two because it is built on the vdW term.

---

## 2. What is stable anyway: the gp5 side

Across **all 28 cluster representatives**, gp5's side of the interface sits in the epitope loop —
**51.1% of total gp5 interface weight falls in B:100–125**, and the per-cluster figure is 42–92%
in every one of the seven clusters. Two residues are contacted in **28/28 models**:

| gp5 residue | models |
|---|---|
| **Arg111** | **28/28** |
| **Gly117** | **28/28** |
| Leu110, Lys118 | 24/28 |
| Gly113 | 23/28 |
| Pro116 | 22/28 |
| Ser112 | 21/28 |
| Leu115, Phe120 | 20/28 |

### This really is the blind run

Worth confirming rather than assuming, because the archive carries no parameter file, and
"every cluster contacts the gp5 loop" is *also* what the one-sided **restrained** run would
produce — there the loop was injected as the active list, and the recovery would be worthless.

Two independent checks say blind:

1. **The active list is not enforced.** Under `job_params_restrained.json` all nine of
   103, 104, 105, 111, 112, 113, 114, 117, 118 carry AIRs, and every model header in this archive
   reports `violations.: 0, 0, 0, …`. A restrained residue with zero violations cannot be absent
   from the interface. Yet:

   | gp5 active-list residue | contacted in |
   |---|---|
   | **Thr105** | **0 / 28 models** |
   | **Asp104** | **2 / 28** |
   | Lys114 | 14 / 28 |
   | Lys103 | 16 / 28 |

   Impossible if these were restrained. 
2. **AIR energies scatter 2.4–84.6 within single clusters** — per-structure random restraints,
   not one shared restraint list.

So the loop identified from escape-mutant genetics was recovered with no information supplied.

### …and it refines the epitope

The blind run recovers **110–118** strongly (R111 and G117 at 28/28, K118 24/28, L110 24/28,
G113 23/28, S112 21/28) and **does not recover 103–105 at all**. The active list was assembled
from escape genetics plus a 40% SASA cutoff; this run says the contact surface is the
**110–118 stretch**, and that 103–105 — solvent-exposed, so it passed the cutoff — is not part of
it. Narrowing the active list to 110–118 is a free improvement to the next restrained run.

**But this is weaker than it looks.** The loop is a **basic patch (net +3, +11.5 per 100 aa) on an
overall acidic protein (gp5 net −13)**, so it is the most electrostatically attractive feature on
gp5's surface, and an elec-dominated blind run should steer there regardless of biology. The
giveaway is that *every* cluster hits the loop while using a *different* aRES face — including the
three clusters (2, 7, 3) that dock it with the RES domain, which no model predicts. Convergence on
one partner's patch while the other partner is smeared is charge steering, not recognition.

---

## 3. What is stable anyway: cluster 6 is a clean bipartite pose

Separating the N-terminal fray (1–5, which carries the charged free N-terminus) from the folded
ribbon is what makes the aRES side interpretable:

| cluster | HS | N-term fray 1–5 | Zn ribbon 6–48 | NTD-1 49–120 | RES 126–447 |
|---|---|---|---|---|---|
| **6** | **−99.9** | **6.2%** | 39.8% | **54.0%** | **0.0%** |
| 5 | −101.3 | 16.2% | 23.2% | 39.4% | 21.2% |
| 1 | −64.9 | **27.0%** | 23.0% | 10.8% | 39.2% |
| 4 | −63.6 | **27.7%** | 20.0% | 12.3% | 40.0% |
| 2, 7, 3 | −88 … −48 | 0% | 0% | 0% | **100%** |

Clusters 1 and 4 looked like N-block poses in aggregate; **~27% of their aRES interface is the
frayed terminus**, so their N-block signal is largely artifact and they should be discounted.
Cluster 6 is the cleanest pose in the run: 94% of its aRES interface is on the genuinely folded
N-block, only 6% on the fray, and **zero contact with the RES domain**.

Decomposing which gp5 surface each aRES module touches in cluster 6:

| aRES module | net charge | gp5 partners | in epitope loop | gp5 side |
|---|---|---|---|---|
| **NTD-1 49–120** | **−16** | R111, S112, G113, L110, K114, G117, P116, K118, L115 | **81%** | 30% K/R, 2% D/E |
| Zn ribbon 6–48 | **+11** (+13 with Zn²⁺) | K359, M273, R288, G360, P284, D358, T357 | 7% | 35% K/R, 10% D/E |

**The acidic module binds the basic loop and the basic module goes elsewhere.** The NTD-1 half is
charge-complementary (acidic −16 against a surface that is 30% K/R and 2% D/E); the ribbon half is
basic-against-basic, i.e. electrostatically unfavourable and probably forced by the random AIRs.

The aRES residues forming the loop-binding face, all in NTD-1:

```
Gln74  Glu89  Tyr92  Glu93  Cys96  Asp97  Asn100 Pro101
His102 Asn103 Asp104 Tyr105 Ile108 Tyr110 Asn111 Asp112
```

and the gp5 loop residues they engage are **107, 108, 110, 111, 112, 113, 114, 115, 116, 117,
118** — which includes **6 of the 9 residues on the independently derived active list**
(111, 112, 113, 114, 117, 118; the three missed are 103, 104, 105) plus 4 passive ones. The
restrained run's active list came from escape genetics and SASA; this run recovered most of it
blind.

### This is the geometry the charge model predicted

The working model recorded in `charge_architecture.tex` — "the basic ribbon binding the DNA and
the acidic NTD-1 binding the polymerase's basic cleft as a DNA mimic" — predicts precisely
cluster 6's arrangement, and it also predicts the ribbon's misbehaviour: **there is no DNA in this
run.** With its intended partner absent, a +13 module has nowhere favourable to go and lands on a
basic gp5 patch at 273–288/357–360 that it should repel.

That the prediction was recorded before this run and that the failure mode is the one the model
implies is worth something. It is still one 8-structure cluster in a non-converged run.

---

## 4. The pre-registered cofold contact was not reproduced — and is dropped

The run's pre-registered readout included one specific residue pair taken from the AF3 cofold: an
acidic residue on the aRES N-block against Ser122 of the gp5 loop. Scored against it:

| | |
|---|---|
| models with the pair in contact (<5 Å) | **0 / 28** |
| closest approach of the pair, in any model | **16.0 Å** (cluster6_2) |

16 Å is not a near miss; it is a different interface. A non-converged run cannot refute a model, so
this is not a refutation — but the run gave the pair **no independent support**, and the pair has
since been **dropped from the project**: not a claim, not a prediction, not a design basis, and not
a crosslink target. More exhaustive docking and MD will nominate contact candidates instead, and
nothing downstream should pre-commit to one.

**This record is kept deliberately.** The readout was registered in advance and it failed, which is
worth more than the original observation ever was — and it is the reason the pair was dropped rather
than re-caveated.

**What survives from this run is at loop level, not residue level:** §2 and §3 above, neither of
which depends on a named pair.

---

## 5. A new and cheaper experiment falls out of this

The persistent salt bridges in the two N-block clusters (<4.0 Å, opposite charge required, counted
per model):

| | | cluster 6 | cluster 5 |
|---|---|---|---|
| aRES **Glu89** – gp5 **Arg111** | NTD-1 | **4/4**, 2.60 Å | — |
| aRES **Glu93** – gp5 **Arg111** | NTD-1 | 2/4, 3.14 Å | — |
| aRES Asp112 – gp5 Lys114 | NTD-1 | 1/4, 2.63 Å | — |
| aRES **Glu93** – gp5 **Lys118** | NTD-1 | — | 2/4, 2.49 Å |
| aRES Asp97 – gp5 Lys118 | NTD-1 | — | 1/4, 2.59 Å |
| aRES Asp112 – gp5 Lys118 | NTD-1 | — | 1/4, 2.68 Å |
| aRES Asp86 – gp5 Arg111 | NTD-1 | — | 1/4, 3.37 Å |
| aRES Asp4 – gp5 Lys359 | *fray* | 4/4, 2.53 Å | 4/4, 2.53 Å |
| aRES Arg29 – gp5 Asp358 | ribbon | — | 4/4, 2.56 Å |

Two aRES-side clusters that agree on almost nothing agree that **gp5 Arg111 and Lys118 are the
basic anchors**, gripped by aRES **Glu89, Glu93, Asp97, Asp112** (and Asp86). Arg111 is contacted
in 28/28 models of the whole run.

That is a specific, independent, plasmid-only experiment that does not depend on docking having
converged, and it targets PF18870 — which has never had a functional test:

- **aRES E89A, E93A** (and the double) through the Fig. 4G co-expression / 4H pull-down assays.
- **gp5 R111A** — predicted to escape, and it is a residue the genetics never nominated.
- Charge reversal **E89K** would be the stronger version if the single alanines are ambiguous.

Note `Asp4 – Lys359` is present in **4/4 models of both clusters** at 2.53 Å. It is the frayed
terminus and should be ignored — it is the most reproducible contact in the run and almost
certainly the least real. A good reminder that reproducibility across models measures the sampling,
not the biology.

---

## 6. Did the zinc "make a big difference"?

**Not assessable from what is in hand.** The pre-Zn blind run's output was never saved to
`results/` — only its parameter file is in the repo — so there is no before/after to compare. What
*can* be said:

- The metal site is **correctly built** (four intended cysteines, 2.08–2.27 Å, tetrahedron
  contracted from 3.67 to 3.53 Å).
- The zinc **makes no direct contact with gp5** anywhere, so any effect is on the aRES fold and on
  the ribbon's +11 → +13 charge, not on the interface itself.
- The `docs/07` §1 caveat — "if docking contacts land on residues 1–48, discount them until the run
  is repeated with Zn modelled" — is now **lifted for residues 6–48**, and replaced by a sharper
  one: discount **1–5**, the charged free terminus, which supplies up to 27% of the aRES interface
  in clusters 1 and 4.
- The run still **did not converge**, so if the pre-Zn run also did not, the honest summary is that
  zinc fixed a real defect in the model without fixing the search.

If the earlier run's archive still exists locally, drop it in and the comparison is one command.

---

## 7. What to do next, in order

1. **Fetch the unclustered HS −154.0 structure** from the server before the two-week deletion. It
   is the best-scoring pose in the run and is absent from the summary archive. If it is an N-block
   pose the reading of this run changes materially; if it is a terminus artifact, that is worth
   knowing too.
2. **AF3 seed replicates of the cofold** (`docs/07` §4), run as a general interface-reproducibility
   check: does *any* NTD-1 contact recur across seeds, and at what ipTM? Let the replicates nominate
   contacts rather than testing a pair chosen in advance (§4).
3. **aRES E89A / E93A + gp5 R111A** through 4G/4H. Cheap, independent of docking convergence, and
   the first functional test aimed at PF18870.
4. **Re-dock with DNA present.** gp5 is a DNA polymerase and the charge model says the ribbon's
   partner is nucleic acid. Its absence is the most likely reason the ribbon half of cluster 6 is
   unphysical. This is also the test that distinguishes "NTD-1 is a protein sensor" from "the
   N-block reads a polymerase·DNA complex".
5. **Replicate seed.** `iniseed: 917` → anything else, same parameters. Cluster 6 holds 8
   structures; whether it reappears is the cheapest check on whether it is real.
6. **Narrow the restrained run's active list to 110–118** before submitting it (§2). The blind run
   recovered that stretch and never touched 103–105, which entered the list on solvent exposure
   rather than on evidence.
7. **Do not** trim 1–5 and re-run as the main fix. Mutating or deleting the terminus changes the
   ribbon's charge. Capping or simply discounting 1–5 in analysis is the right treatment.

### Registered in advance, for the replicate run

- Cluster 6 reappears with NTD-1 dominant and RES near zero → the pose is real, proceed to 3 and 4.
- An N-block cluster reappears but dominated by the ribbon with NTD-1 minor → the DNA-mimic model
  is wrong in its assignment, and the ribbon is the sensor.
- Nothing reproduces → the blind approach is exhausted at this system size. Stop running it, and
  put the compute on the restrained run and the DNA-containing complex.
