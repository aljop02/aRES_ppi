# Blind HADDOCK runs: results, and the P1 verdict

Analysis of the two finished HADDOCK 2.4 runs (jobs **779823** and **780456**), from the uploaded
plot archives and cluster-summary downloads.

**P1 was registered in advance** (`README.md`: *"Commit to scoring the blind run against P1 before
looking at the output"*). It is scored here as registered. The registered answer turns out to be
**"no conclusion"** — not because P1 failed a test, but because neither run converged well enough to
constitute one. §7 holds the pre-registered decision rule and why this is the outcome it names.

---

## 0. The short version

| | |
|---|---|
| **P1 (NTD-1 binds the gp5 epitope loop)** | **Not supported by either run.** No cluster places aRES NTD-1 against gp5 102–122. |
| **Did either run converge?** | **No.** ~94% of water-refined models are unclustered in both; the largest cluster is 3 models. |
| **Is the top-scoring pose trustworthy?** | **No.** Its score is 63–70% electrostatics, and the Zn(Cys)₄ ribbon carries no zinc. |
| **Is this evidence *against* the hypothesis?** | **No.** A non-converged ab-initio run is uninformative, not falsifying. See §7. |
| **The restrained run** | **`2_ares_T7` was missing, then ran and failed.** Job 780457 produced no complexes at all — the partners never touched. §6, §10. |

The honest summary: these two runs do not discriminate between hypotheses. They should be reported
as a non-converged ab-initio result and not interpreted pose-by-pose.

---

## 1. Which run is which

Both plot archives lacked run labels. They are identified by matching the cluster numbering in the
summary downloads and the analysed-model count against the recorded parameter files.

| | **job 779823** | **job 780456** |
|---|---|---|
| Server run name | `aRES-T7` | `1_aRES_T7` |
| Plot archive | `plots1.tgz` | `plots.tgz` |
| Parameter file in repo | — (superseded v1) | [`job_params_blind_v2.json`](../analysis/docking/blind_run/job_params_blind_v2.json) ✅ |
| `ranair` | `true` (inferred, §4) | `true` (recorded) |
| Water-refined models analysed | **200** | **400** |
| Clusters found | 5 | 12 (top 10 downloadable) |
| Cluster models shipped | 11 | 20 |
| Best HADDOCK score | **−165.03** (model 157) | **−168.27** (model 108) |
| aRES construct (chain A) | **1–459** — tail present | **1–447** — tail trimmed |
| Date of coordinates | 13-Oct-2026 | 13-Oct-2026 |

`job_params_blind_v2.json` declares `runname: 1_aRES_T7`, `ranair: true`, `anastruc_1: 400` — all
three match job 780456, so **780456 is the widened-funnel blind run** the repo documents.

Job 779823 matches neither recorded parameter file and analysed 200 models, so it is the earlier
blind run (`ares_T7_a`, reviewed in `docs/07` from the since-superseded `job_params.json`). Note it
analysed **200** models, not the 50 described in
[`blind_run/README.md`](../analysis/docking/blind_run/README.md) — so the run that actually
executed was not the 50-structure configuration that README criticises.

The aRES construct settles the ordering independently: 779823 carries the **full 1–459 chain**,
while 780456 stops at **447**. `blind_run/README.md` prescribed "bump sampling to 10000, trim the
tail" for v2, and 780456 is the run with the tail trimmed. The two runs therefore differ in *two*
variables, not one — sampling depth **and** construct — so they are not a clean pair for judging
whether the wider funnel helped.

**In both complexes chain A = aRES and chain B = gp5 (1–704).** This matches the usage in both
run READMEs but is the *opposite* of the example in the top-level `README.md`
(`--domain B:...:NTD1`), which has been corrected.

---

## 2. The P1 verdict: not supported

P1 as registered: **the aRES N-terminal domain (NTD-1, 49–120) contacts the gp5 epitope loop
(insertion1, 102–122).** Both halves must hold.

Interface composition aggregated over all shipped cluster models, from
[`analysis/scripts/interface_contacts.py`](../analysis/scripts/interface_contacts.py):

| aRES side (chain A) | 779823 | 780456 |
|---|---|---|
| Zn ribbon (1–48) | 15.0% | 9.7% |
| **NTD-1 (49–120)** | **32.7%** | **19.0%** |
| hinge (121–125) | 0.0% | 0.0% |
| RES (126–447) | **45.6%** | **71.3%** |
| tail (448–459) | 6.8% | — (trimmed) |

| gp5 side (chain B) | 779823 | 780456 |
|---|---|---|
| **insertion1 (102–122)** | **8.6%** | **5.4%** |
| TBD (258–333) | 33.7% | 17.2% |
| 410–704 (coarse bin, C-terminal half) | 28.9% | 47.3% |

The aRES side of the interface is **majority RES domain in both runs**, and the gp5 epitope loop
accounts for 5–9% of interface weight — roughly its share of the gp5 surface, i.e. no enrichment.

In run 779823 the untrimmed tail is itself in the interface (F453, D455, L456, D458, F459, each in
2/11 models). A disordered C-terminal tail acting as an interface is a sampling artifact, and it is
the reason the tail was trimmed for 780456 — so 779823's interface statistics are the less reliable
of the two even before §3 and §4.

### Per cluster, no cluster satisfies both halves

Each cluster is a single, essentially pure patch (full table in
[`results/haddock_blind_clusters.tsv`](../results/haddock_blind_clusters.tsv)):

| run | cluster | best HS | aRES side | gp5 side | P1? |
|---|---|---|---|---|---|
| 779823 | 1 | −147.73 | NTD-1 68.6% (+Zn 31.4%) | **TBD 88.7%**, ins1 0% | ✗ wrong gp5 surface |
| 779823 | 5 | −104.69 | RES 100% | ins1 30.2% | ✗ wrong aRES domain |
| 780456 | 3 | −67.93 | NTD-1 100% | ins1 0%, 410–704 94.4% | ✗ wrong gp5 surface |
| 780456 | 6 | −51.23 | NTD-1 100% | ins1 0%, 410–704 100% | ✗ wrong gp5 surface |
| 780456 | 10 | −116.18 | RES 100% | ins1 25.0% | ✗ wrong aRES domain |
| 780456 | 11 | −100.02 | **Zn ribbon 90.6%** | TBD 78.4% | ✗ discounted, §4 |

The clusters that get the aRES domain right get the gp5 surface wrong, and vice versa. **No cluster
in either run puts NTD-1 against insertion1.**

Cluster 11 of run 780456 is 90.6% Zn ribbon, which falls squarely under the rule pre-registered in
[`docs/07` §1](07_docking_run_review.md) ("Why it matters", point 3): *"If docking contacts land on
residues 1–48, discount them until the run is repeated with Zn modelled."* It is discounted.

---

## 3. Neither run converged

This is why §2 should not be read as evidence against the hypothesis.

| | 779823 | 780456 |
|---|---|---|
| Models not assigned to any cluster | **189 / 200 (94.5%)** | **376 / 400 (94.0%)** |
| Largest cluster | **3 models** | **2 models** |
| Best-scoring model | unclustered singleton | unclustered singleton |
| Models within 10 Å i-RMSD of the best model | **1 / 205** | **2 / 410** |
| Max residue–residue contact frequency | 0.27 (3/11) | **0.10 (2/20)** |

Denominators differ by design: **200 / 400** are distinct water-refined models, while **205 / 410**
are rows in the analysis plots, which carry one extra per-cluster mean row per cluster.

No reference structure was supplied, so HADDOCK measured i-RMSD/FCC **against the best-scoring
model**. Only that model itself sits at i-RMSD 0 / FCC 1.0; the next-nearest pose is 19 Å away in
779823 and 25 Å in 780456. There is no funnel — the score–RMSD plots show a flat band, not a basin.

**The wider funnel did not help.** Doubling the analysed set (200→400) produced more and smaller
clusters (5→12) and dropped the best reproducible contact frequency from 0.27 to 0.10. That is the
signature of sampling a featureless energy landscape: more sampling finds more distinct mediocre
poses, not more copies of one good pose.

Strictly, the two runs differ in both sampling depth and aRES construct (§1), so this is not a
controlled comparison of funnel width — but neither run converged, which is the point that matters.

⚠️ **Open discrepancy.** The shipped clusters contain 2–3 models each. HADDOCK ships up to the top
4 models per cluster, so these clusters genuinely are that small — yet both recorded parameter files
set `clust_size: 4`, which should have discarded them. Worth checking against the server's run page
before the results expire.

---

## 4. The score is electrostatics, and there is no zinc

The HADDOCK water-stage score reconstructs exactly from the reported terms
(max |error| 5 × 10⁻⁵ over all 615 plotted models):

```
HADDOCK score = Evdw + 0.2 * Eelec + Edesolv + 0.1 * Eair
```

**Eair is exactly 0 for every model in both runs** — in the plots ("No Eair was calculated for this
run") and in every PDB `REMARK energies` line. That is the `ranair: true` signature: random AIRs are
dropped before water refinement, so nothing restraint-derived survives into the final score. It also
confirms job 779823 was ab-initio, and it is the diagnostic that proves the restrained run is absent
(§6).

With Eair gone, the ranking is left to three terms, and electrostatics dominates:

| | 779823 | 780456 |
|---|---|---|
| Share of top model's score from 0.2·Eelec | **62.8%** | **69.6%** |
| Median share across all models | 43.5% | **55.1%** |
| Most negative Eelec | **−517.9** | **−585.4** kcal/mol |
| Models with Eelec < −400 | 5 / 205 | 9 / 410 |
| Models with *unfavourable* (positive) Edesolv | 33 / 205 (16%) | **192 / 410 (47%)** |
| Median BSA of cluster models | 892 Å² | 1254 Å² |

Three things are wrong here at once:

1. **Eelec of −500 to −585 kcal/mol is not a protein–protein interface**, it is a charge-pairing
   artifact. Cluster 12 of run 780456 is the extreme case: Eelec −497 with Edesolv **+18.9** — a
   large electrostatic reward paid for by burying charged groups that would rather stay solvated.
2. **Buried surface area is too small.** Median 892 Å² / 1254 Å², and 15 of 31 cluster models are
   below 1200 Å². HADDOCK reports BSA summed over both partners, so 892 Å² is ~450 Å² per side —
   crystal-contact territory. An average specific PPI buries roughly 1500–2000 Å² total.
3. **The Zn(Cys)₄ ribbon has no metal.** No `HETATM` record appears in any model. `docs/07` found
   C8/C11/C31/C34 in near-ideal tetrahedral geometry; without Zn²⁺ that site is four uncoordinated
   thiolates plus a +2 charge hole, which is exactly what an electrostatics-driven search will home
   in on.

**[INTERPRETATION]** gp5 is a DNA polymerase, so it carries a strongly basic DNA-binding cleft; aRES
presents acidic patches and an unneutralised Zn site. Blind docking with no restraints, no counter-
ions and no metal will preferentially pair those, which accounts for both the huge Eelec values and
the absence of a funnel. The concrete test: overlay the top poses on the DNA path in **1T7P** and
check whether the aRES contact footprint on gp5 coincides with the nucleic-acid groove. If it does,
these poses are charge artifacts and nothing more.

---

## 5. The one feature both runs share

aRES residues **138–144** (M138, K139, E140, T141, R142, Y143, F144) are the top-ranked interface
patch in *both* runs — frequency 0.55 (6/11) in 779823 and 0.40 (8/20) in 780456 — despite the two
runs differing in construct and sampling depth and sharing not one cluster between them.

It is a mixed-charge, aromatic-containing patch immediately C-terminal to the hinge, at the start of
the RES block. Two readings, not distinguishable from these data:

- it is the most exposed sticky patch on the isolated aRES model, and will attract any partner in
  any charge-driven docking run — the null explanation, and the one to beat;
- it is a genuine hot spot.

Treat it as a flagged observation only. It is **not** in NTD-1, so it is not relevant to P1 either
way. If it survives a run with Zn modelled and real restraints, it becomes worth a mutant.

---

## 6. ⚠️ The restrained run has not been run

[`job_params_restrained.json`](../analysis/docking/onesided_run/job_params_restrained.json) declares
`runname: 2_ares_T7` with `ranair: false` and 9 active gp5 residues. **Neither upload is that run:**

- the two server run names are `aRES-T7` and `1_aRES_T7`, neither of which is `2_ares_T7`;
- `ranair: false` with real AIRs would keep a restraint term through water refinement, so **Eair
  would be non-zero**. It is exactly 0 in all 615 models.

So the one-sided run — the one designed to answer P1 by pinning the gp5 epitope and leaving the aRES
side free — is still outstanding. **That run is the experiment; these two were the control that shows
blind docking cannot answer the question.**

When it comes back, **check Eair first.** If the `air_boxplot.html` again says *"No Eair was
calculated for this run"*, the restraints did not take effect and the run is another blind run
regardless of what was submitted.

> **Update, same day.** It was submitted as job **780457** and came back having never docked
> anything. The Eair check above fired immediately and was the right diagnostic. See §10.

---

## 7. What these runs do and do not license

The applicable decision rule was written down in advance, in
[`docs/07` §3](07_docking_run_review.md):

> - Blind converges on gp5 ~100–125 **and** aRES NTD-1 → the strongest computational result
>   available in this project. Nothing went in; it came out anyway.
> - Blind converges somewhere else, reproducibly across seeds → a real alternative hypothesis worth
>   taking seriously.
> - **Blind does not converge → no conclusion.** Not a strike against the hypothesis. Move on.

**The outcome is branch three.** §3 establishes non-convergence on every available measure, so the
pre-registered conclusion is *no conclusion* — and that was committed to before the output was seen.

Note that the `≥60% of aRES interface weight in 1–120` threshold in
[`onesided_run/README.md`](../analysis/docking/onesided_run/README.md) is the readout for the
**restrained** run, not these. For the record the blind runs give 47.7% (779823) and 28.7% (780456)
in 1–120, and aRES **D104 is essentially absent** from the interface — 0/11 models in 779823 and
1/20 in 780456. Those numbers are reported for completeness, not applied as a test: that table
assumes a run that localised something.

**Do not say:** "docking does not support NTD-1 as the sensor." Both runs failed at the sampling
stage, before reaching the point where they could test anything. §2 is an absence of signal in a run
with no signal anywhere.

**Do say:** blind ab-initio docking of aRES against gp5 alone does not converge, its scoring is
dominated by an electrostatic artifact traceable to an unmodelled Zn site, and the question needs
restraints — which is what the one-sided run was designed to supply.

**The cofold result stands untouched.** The AF3 aRES×gp5 cofold contact (gp5 Ser122 OG ··· aRES
Asp104, 3.06 Å; `docs/05` §2) is independent of these runs and is not weakened by them.

---

## 8. Next steps, in order

1. **Download the two best-scoring models before they expire.** Both are unclustered singletons
   (779823 model 157, HS −165.03; 780456 model 108, HS −168.27) and are therefore **not** in the
   cluster-summary archives. HADDOCK deletes server results after two weeks.
2. **Resubmit `2_ares_T7` with the explicit `ambig.tbl`** — the first attempt (job 780457) failed
   with an empty restraint table and never docked anything (§10). Run the 100-structure test job
   first and verify Eair ≠ 0 and BSA > 0 before committing the full sampling.
3. **Model the Zn²⁺** into the Zn(Cys)₄ site, or exclude 1–48 from the interface, before any further
   blind docking. Until then, every blind result inherits the artifact in §4.
4. **Dock gp5·TrxA** (T1.3 in `docs/03`): take 1T7P, strip the DNA, keep thioredoxin. The TBD hits
   in §2 are most likely an artifact of a surface that is covered in the activating species — but
   the diagnostic table in [`docs/02` §4](02_hypothesis_NTD1_sensor.md) cannot be read until the run
   with TrxA present exists, because gp5-alone docking cannot distinguish "aRES binds the TBD" from
   "the TBD is artificially exposed."
5. **Overlay the top poses on the 1T7P DNA path** (§4) — a one-afternoon check that would settle
   whether the blind poses are charge artifacts.

None of this changes the bench priorities in `docs/05` §4 or `docs/08`. The decisive experiments
there (gp5 S122A, aRES D104A, the λ P/O co-expression, the UC5↔GN04826 chimeras) do not depend on
any docking result.

---

## 9. Reproduce

Derived tables, kept in `results/` per the `.gitignore` policy of not committing raw run output:

| File | Contents |
|---|---|
| [`results/haddock_blind_clusters.tsv`](../results/haddock_blind_clusters.tsv) | 15 clusters × energies, BSA, and the full domain/ROI breakdown |
| [`results/haddock_blind_contacts_779823.tsv`](../results/haddock_blind_contacts_779823.tsv) | per-residue interface contact frequency, 11 models |
| [`results/haddock_blind_contacts_780456.tsv`](../results/haddock_blind_contacts_780456.tsv) | per-residue interface contact frequency, 20 models |

From the unpacked cluster-summary download:

```bash
python3 analysis/scripts/interface_contacts.py <run_dir>/*.pdb --chains A B \
    --domain A:1-48:Zn_ribbon --domain A:49-120:NTD1 \
    --domain A:121-125:hinge --domain A:126-447:RES --domain A:448-459:tail \
    --roi B:102-122:gp5_insertion1 --roi B:258-333:gp5_TBD --roi B:410-704:gp5_polC \
    --tsv results/haddock_blind_contacts_<job>.tsv
```

Keep `--domain A:448-459:tail` even for the trimmed construct: it is what makes the breakdown sum
to 100%, and its absence is how the untrimmed run was spotted.

Per-model energies and BSA come from the `REMARK energies` / `REMARK buried surface area` lines of
the cluster PDBs; the plotted score/RMSD tables are embedded as Bokeh JSON in the `*.html` plot
files.

---

## 10. Job 780457 (`2_ares_T7`) — the restrained run failed to dock anything

Submitted after §6 was written, returned the same day. It is **not a docking result**: the two
partners never came into contact, so nothing in it bears on aRES, gp5, or P1.

### What came back

| | |
|---|---|
| Structures analysed | 400 |
| Clusters | **400, every one of size 1** |
| Van der Waals energy | **exactly 0.0** |
| Electrostatic energy | **exactly 0.0** |
| Restraints violation energy (Eair) | **0.0** |
| Buried surface area | **−999999** (HADDOCK's error sentinel) |
| HADDOCK score | equals the desolvation energy exactly (−6.7 = −6.7, −5.9 = −5.9) |

### The measurement that settles it

Chain positions in the ten downloadable models:

| model | COM chain A (aRES) | COM chain B (gp5) | separation |
|---|---|---|---|
| cluster391_1 | (−167.0, 0.2, −0.2) | (108.0, −0.3, 0.1) | 275.0 Å |
| cluster393_1 | (164.6, −1.2, 0.1) | (−106.5, 0.7, 0.0) | 271.1 Å |
| cluster399_1 | (−162.0, 2.1, −0.4) | (104.7, −1.4, 0.4) | 266.7 Å |
| … all ten | | | **266–283 Å** |

Both chains are present and complete (aRES 1–447, gp5 1–704). They sit on the x-axis with y,z ≈ 0 —
**their untouched starting positions.** Only the sign of x varies, which is `rotate180_0: true`.
Minimum inter-chain atom distance in cluster391_1 is **217 Å**.

### Root cause: the AIR table was empty

With `ranair: false`, `cmrest: false` and `surfrest: false`, the restraint energy is the *only* term
that can bring two separated molecules together during it0 rigid-body minimisation. At 217 Å, vdW
and electrostatics are far beyond their cutoffs, so there is not even a gradient pointing the right
way. No restraints → no association → the molecules stay exactly where they were placed.

The restraint **definition** in `job_params_restrained.json` is sound, which is what makes this
worth recording rather than just rerunning:

| check | result |
|---|---|
| gp5 active residues 103–118 present in the uploaded PDB | ✅ all 9 |
| aRES passive residues present in the uploaded PDB | ✅ all 304 |
| segid / chain consistency | ✅ `A`/`B` on both |
| `inter_mat` | ✅ all 1.0 — the molecules are allowed to interact |
| scoring weights | ✅ standard |
| `tblfile` / `unambigtblfile` | **`None`** — no uploaded restraint file |
| `partners.1.activereslist` | **`[]`** — zero active residues on aRES |

So the server was asked to synthesise the AIRs from the active/passive lists and produced none.
HADDOCK builds restraints *outward from active residues*; with zero on aRES, only gp5's 9 active
residues should have generated any, and none appeared. The "Submit file → Full parameter file"
route may also bypass the restraint builder altogether.

**[INTERPRETATION]** The zero-active-residues configuration on partner 1 is the likely trigger. The
one-sided design itself is standard and correct — active epitope on the known side, solvent-
accessible surface as passive on the unknown side — it just cannot survive being handed to a
generator that wants active residues on both partners.

### Fix

Restraints are now written explicitly rather than inferred, by
[`analysis/scripts/make_ambig_tbl.py`](../analysis/scripts/make_ambig_tbl.py), and committed as
[`analysis/docking/onesided_run/ambig.tbl`](../analysis/docking/onesided_run/ambig.tbl) — 9
restraints, one per gp5 active residue, each targeting all 304 aRES surface residues:

```bash
python3 analysis/scripts/make_ambig_tbl.py \
    --from-params analysis/docking/onesided_run/job_params_restrained.json \
    --source-segid B -o analysis/docking/onesided_run/ambig.tbl
```

Pre-flight checks before the next submission are in
[`onesided_run/README.md`](../analysis/docking/onesided_run/README.md): confirm the table is
non-empty, verify segids survive upload, set `noecv: false` (50% random restraint removal is
reckless with only 9 AIRs), and burn ten minutes on a 100-structure test job to prove `Eair ≠ 0`
and `BSA > 0` before committing 5000 structures.

### The transferable lesson

`Eair = 0` is the diagnostic that separates **"docked badly"** from **"never docked"**, and it has
now caught two different failures: random-AIR dropout in the blind runs (§4) and a missing restraint
table here. Check it before reading anything else off a HADDOCK run. A results page full of
plausible-looking cluster numbers is not evidence that docking occurred — these 400 clusters had
Z-scores, RMSDs and rankings, and not one of them was a complex.
