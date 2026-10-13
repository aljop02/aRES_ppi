# Structural characterisation at EMBL Hamburg — rotation plan

**How this was produced.** Five technique experts (crystallography, cryo-EM, biophysics, protein
production, interface mapping) planned in parallel, each proposal was stress-tested by an adversarial
feasibility critic, and the surviving experiments were merged into the plan below. Facility details
came from web search only — `embl.org`, `embl-hamburg.de` and `cssb-hamburg.de` are all blocked by
this session's egress policy — so **everything about instruments and access is unverified and
collected in §12 for you to confirm on arrival.**

## Independently verified against the repository's own coordinates

I re-ran the structural claims against `analysis/docking/*/job_params_*.json` (partner 1 = aRES
chain A, 1–447; partner 2 = gp5 chain B, 1–704). All of them hold:

| Claim | Status |
|---|---|
| Zn ribbon C8, C11, C31, C34 | ✅ all Cys |
| **Two further cysteines in the N-block: C78, C96** | ✅ six Cys in 1–120, not four |
| **No lysine anywhere in 49–120** | ✅ all ten N-block lysines are in the ribbon (K3, K5, K10, K15, K23, K35, K37, K39, K44, K48) |
| Five methionines in 1–120 (M1, M18, M49, M50, M60) | ✅ SeMet is viable |
| D104 is Asp | ✅ |
| R182, R184, R215, E236 | ✅ all match the paper |
| **Residue 234 is Ile, 293 is Ser** | ⚠️ the paper's defence-null list says **S234A** and **V293A** |

### Two corrections to earlier documents in this repo — both mine

**1. The RES Pfam accession is `PF08808`, and my earlier "correction" of your note was wrong.**
I told you `PF08808` was a typo for `PF08808`. It is not. The paper itself is internally
inconsistent: **Fig. 1B labels the domain `PF08808 RES`**, while the STAR Methods text repeatedly
says `PF08808`. `PF08808` is a histone-like family (CBFD_NFYB_HMF), not RES. Your original note and
your PAE figure had it right. This matters concretely: an HMM search built on `PF08808` returns the
wrong family, which would silently corrupt the homolog collection for the conservation analysis in
`docs/03` T4. Confirm on InterPro when the network allows — I could not reach it from here.

**2. The 234/293 mismatch looks like a transposition, not a different isolate.**
The model has **Ile234** and **Ser293**; the paper lists **S234A** and **V293A**. Note the pattern:
the model's Ser is at 293 and its 234 is a branched aliphatic (Ile, conservative to Val). That is
what a swapped pair of residue numbers looks like — i.e. the paper may mean S293A and V234A. The
alternative is that the model is a different UC5 isolate. **This is upstream of every mutant in the
plan and of the UC5↔S_34AVR map**, so resolve it against the primary sequence record before
ordering anything. Four of six positions verify, so this is a localised problem, not a wholesale
frame error.

## ⚠️ D104 is NOT load-bearing — read this before §2

The plan below was written treating the AlphaFold3 cofold contact (gp5 Ser122 → aRES Asp104) as a
working model. **That is too strong.** It is a single observation, in one seed, at poor ipTM. It has
not earned the weight several sections give it, and no resource should be committed to it before
independent evidence arrives. The docking runs in `analysis/docking/` are exactly that evidence
being gathered — do not pre-empt them.

**What dies if the cofold is wrong:** D104 as the specific partner; the pre-registered crosslink
target set (D104, D106, D97, E98, D112); aRES D104A/D104N as constructs; the "D106 scaffold /
S122 contact" two-step model in `docs/05` §2, which was derived from the cofold's own geometry.

**What is untouched, because it never came from the cofold:**

| | Source of evidence |
|---|---|
| The gp5 epitope is the ~100–125 loop | **Phage genetics** — two independently isolated escape mutants (S122P, D106G), each escaping *both* the UC5 and S_34AVR aRES. The hardest evidence in the project |
| NTD-1 is the sensor module | The authors' own statement + domain architecture (PF18870 fused to RES, HEPN, ABC ATPase) |
| The Zn(Cys)₄ ribbon | Monomer model geometry, pLDDT 90.8 |
| Domain boundaries 1–120 / hinge / 126–447 | Monomer PAE |
| NTD-1 has no lysine, so BS3/DSS cannot report it | A fact about the sequence |
| GN04826 as the natural non-binder | Fig. 1B |
| **All of Track D** | One protein; no interaction required |

**And the experiments that matter most were never D104-dependent:** does the isolated N-block bind
gp5·TrxA at all; does it compete in trans; does the gp5 98–132 peptide activate; does
gp5(Δ118-123)·TrxA still activate; what is the Kd. Each tests "is NTD-1 the sensor" and "is the loop
the epitope" **without naming a residue pair.** Those stay exactly as planned.

**Consequences applied below:** P2 is restated at the level the evidence supports; aRES D104A/D104N
leave the week-1 gene order; XL-MS becomes a *discovery* experiment across the whole acidic face of
the N-block rather than confirmation of a pre-registered pair (the prediction is still written down
first — it is just no longer the design basis); and Decision Point 1 asks whether *any* NTD-1
contact reproduces, not whether that one does.

---

## The plan

*All residue numbers are **E. coli UC5** unless marked S_34AVR. Protein must be made from the K. pneumoniae S_34AVR homolog — UC5 aRES is not solubly expressible — so a numbering map is a week-1 prerequisite, not a formality.*

---

## 0. Bottom line

A rotation cannot deliver a structure of the 143 kDa aRES·gp5·TrxA complex. It can plausibly deliver **one of two publishable things**, and the plan is built so that both are live until a decision point kills one:

| | Deliverable | Realistic odds in 12–14 wk |
|---|---|---|
| **Track D** | First experimental structure of any PF18870 domain (apo N-block, ~14 kDa), by MR on the AlphaFold model with Zn anomalous validation | **~30–35%** |
| **Track I** | A measured Kd + a residue-level interface restraint set (zero-length XL-MS) + the escape-mutant and natural-non-binder controls | **~40–45%** |
| Either | At least one of the two | **~60%** |

The two tracks share one protein prep and then diverge. **Track D does not depend on the interaction work succeeding at all** — that is its main virtue, because the interaction hypothesis rests on a single AlphaFold3 contact.

Three things will actually decide the rotation, and none of them is a beamline:
1. **Does the isolated N-block express solubly, fold, and load zinc?** (gates all of Track D)
2. **What is the Kd?** (sets which technique in Track I is even applicable — every method has a concentration window and nobody has measured it)
3. **Does recombinant gp5·TrxA exist at ≥0.5 mg by week 6?** (gates all of Track I; the commercial NEB enzyme cannot supply it)

---

## 1. Two goals, and the asymmetry between them

**The interaction (aRES ↔ gp5·TrxA).** High value, high dependency. Needs two proteins, one of them an 80 kDa phage polymerase co-expressed with TrxA, a measured Kd, and an interface chemistry that can physically report the predicted contact. Four independent gates.

**The domain (NTD-1 / PF18870 itself).** Lower ceiling per unit novelty — but **there is no experimental structure of any PF18870 domain anywhere**, so an apo N-block structure, or even a rigorous solution characterisation of it, is a first. One protein, one gate (does the fragment behave?).

**The fallback relation is directional:** if the interaction work collapses (gp5 won't express, Kd is >50 µM, the N-block doesn't bind), Track D continues untouched. The reverse is not true — if the N-block is insoluble, Track D dies and Track I continues on full-length aRES. So **the N-block expression screen and the gp5·TrxA co-expression must start in the same week**, in parallel, and neither may be made to wait on the other. Serialising them is the single most likely way this rotation ends with protein in a freezer and no data.

---

## 2. The five predictions the whole project rests on

Every one of these is currently a prediction. For each: what kills it, and how cheaply.

| # | Prediction | Current evidence | Experiment that kills it | Cost |
|---|---|---|---|---|
| **P1** | The AlphaFold3 cofold interface is real | One cofold, one pose | **Interface-confidence audit**: ipTM, interface PAE, and per-seed frequency of the D104–S122 contact across all seeds. If it appears in 1 of 5 seeds, stop and rethink before ordering a primer | **Free, day 1** |
| **P2** | The gp5 100–125 loop contacts the aRES N-block *(restated — the earlier form, "Ser122 contacts Asp104", overstates one low-ipTM seed)* | Escape-mutant genetics for the loop; the residue pair rests on one cofold only | **Zero-length XL-MS** (EDC/DMTMM): gp5 K103/K114/K118 → an aRES acidic residue in 49–120. Plus gp5(Δ118-123)·TrxA activation test — a *published, catalytically competent* deletion construct (1T7P/1TK0 carry DEL(118-123)) that removes G121/S122 outright. If Δ118-123 still activates aRES, S122 is not the contact | XL-MS ~100 µg; Δ-construct is one clone |
| **P3** | NTD-1 is the sensor module (not the RES block, not the Zn ribbon) | Domain architecture + the authors' own statement; *not* the cofold | **Does the isolated N-block bind gp5·TrxA?** MST/nanoDSF-FoldAffinity + protease protection. And: does excess N-block *compete* in trans against gp5-driven activation of full-length aRES, on the established NAD/NADH-Glo assay? | ~1 µg labelled N-block per titration; one plate |
| **P4** | aRES carries a structural Zn(Cys)4 site | AlphaFold SG–SG geometry only (3.55–3.81 Å, ideal 3.76; no metal in the model) | **ICP-MS/OES stoichiometry + native ESI-MS (+63.4 Da) + Zn K-edge XANES at PETRA III P64.** XAS gives coordination number, ligand type (S vs N/O) and Zn–S distance directly | ICP ~0.2 mg; XAS is on-campus |
| **P5** | aRES is a monomer and the complex is 1:1 | Assumed throughout | **SEC-MALS or AUC on aRES alone, before any stoichiometry is quoted.** The only solved RES-NADase structure (6GW6, P. putida RES–Xre, 2.2 Å) is built on a central **RES toxin dimer**. If aRES dimerises, apo aRES is ~100 kDa — unresolvable from 92 kDa gp5·TrxA by mass photometry — and every mass, mask volume and particle-count figure in the cryo-EM plan is wrong | One SEC run |

**P1 and P5 are free and must be done in week 1.** They are the two cheapest ways to discover that a downstream month was going to be wasted.

---

## 3. What I verified from the repository's own coordinates

Both monomer models are present, embedded as `partners['1']['raw_pdb']` (aRES, **chain A**, residues 1–447) and `partners['2']['raw_pdb']` (gp5, **chain B**, 1–704) inside `analysis/docking/onesided_run/job_params_restrained.json` and `.../blind_run/job_params_blind_v2.json`. Note the chain convention: **A = aRES, B = gp5** (confirmed against `analysis/scripts/interface_contacts.py` usage and both docking READMEs). Any script written the other way round silently reports the wrong chain.

Four results that change the plan:

**(a) NTD-1 is a lysine desert — this kills BS3/DSS as an interface probe.**

Residues 49–120 contain **no lysine at all**. The only basic residues are R69, R109, R117. All ten N-block lysines are in the Zn ribbon: K3, K5, K10, K15, K23, K35, K37, K39, K44, K48.

Nearest Lys of any kind to Asp104 (Cα–Cα): K48 **21.3 Å**, K44 23.3, K3 24.6, K5 28.0 — all ribbon — then RES-block K274 **29.3 Å**, K281 **32.2 Å**.

Consequences: (i) a Lys–Lys crosslinker can *never* place the gp5 loop "on NTD-1"; (ii) **a BS3 crosslink from the gp5 loop to K274 or K281 sits inside the ~30 Å convention and therefore does not refute the NTD-1 hypothesis** — do not read it as refutation; (iii) the gp5 side is fine (K103 at 8.2 Å from S122, K118 at 8.3, K114 at 13.6), so the problem is entirely on the aRES side.

**(b) There is a genuine acidic patch around Asp104**, which is the real reason zero-length chemistry suits this interface: D106 (5.6 Å), D97 (13.2), E98 (13.6), D112 (16.8), E93 (17.7), E62 (17.0). Pre-register these as the EDC/DMTMM target set.

**(c) SeMet is viable — this was an open question and the answer is good.** The N-block 1–120 contains **five methionines (M1, M18, M49, M50, M60)**, i.e. one per 24 residues. Excellent anomalous density for a 14 kDa fragment. 49–136 retains M49/M50/M60. So experimental phasing has a second route if zinc turns out substoichiometric.

**(d) There are two non-ribbon cysteines in the N-block: C78 and C96.** Six cysteines total in 1–120, not four. This raises oxidation and disulfide-scrambling risk above what the Zn-ribbon story implies, it confounds any thiol-directed zinc-release chemistry (PMPS/MMTS), and it means the 49–136 "no-ribbon" construct still carries two free cysteines with no metal to coordinate them. Keep TCEP at **1–2 mM**, not 0.5 mM, and degas buffers. (Also: W2 and W85 are in 1–120, so nanoDSF intrinsic fluorescence works on the N-block.)

**(e) A numbering problem inside the UC5 frame itself.** In this model, residue 234 is **Ile** and residue 293 is **Ser**, whereas the project's defence-null list reads **S234A** and **V293A**. R182 (Arg), R184 (Arg), R215 (Arg), E236 (Glu) and C8/C11/C31/C34 (all Cys) check out. So four of six catalytic positions verify and two do not. Either the model is of a different UC5 isolate or two mutants were transcribed in a shifted frame. **The UC5↔S_34AVR map cannot be built until the UC5 frame is reconciled against the primary sequence record.** This is a week-1 item and it is upstream of every mutant in the plan.

**(f) The architecture is extended, which makes one SAXS experiment sharp.** All-atom Rg of aRES 1–447 is **31.9 Å**; N-block (1–120) Rg 15.5 Å, RES block (126–447) Rg 24.6 Å, **centroid separation 50.5 Å**. A compact 52 kDa globule would be ~25 Å. The N-block is **27% of the mass of aRES alone** but only ~10% of the 143 kDa complex — which is exactly why SAXS is well-posed on apo aRES and ill-posed on the complex (§8.5).

---

## 4. The gate structure

```
WEEK 0-1  DRY + ORDERS + EMAILS  (nothing is gated on anything)
  ├─ P1 cofold interface audit  ──┐
  ├─ P5 aRES oligomeric state   ──┤── if either fails, re-plan before spending money
  ├─ UC5 frame reconciliation + S_34AVR map
  ├─ Gene orders out (aRES boundary series + gp5+trxA operon)
  └─ Emails: Sorek lab plasmids / Chiu operon / SPC iLAB / P64 / P14 / cryo-EM

                    ↓ (2-3 wk synthesis clock runs in background)

GATE 1  EXPRESSION      small-scale parallel screen, lysate SDS-PAGE + 96-well IMAC
  ├─ N-block soluble? ───────────────→ TRACK D opens
  └─ gp5·TrxA ≥0.5 mg? ─────────────→ TRACK I opens
     (these are independent; one can fail without the other)

GATE 2  FOLD           nanoDSF (cooperative Tm) + far-UV CD + DLS + ICP-MS Zn
  │   a cooperative Tm alone is NOT a fold gate for a novel fold — CD is mandatory
  └─ fails → boundary redesign via limited proteolysis; do not spend protein on plates

GATE 3  Kd            nanoDSF FoldAffinity → MST → (ITC only if affordable)
  ├─ Kd < 1 µM   → full Track I works, incl. SEC-SAXS and crosslinking at low conc.
  ├─ Kd 1-20 µM  → XL-MS and protease protection survive; mass photometry and
  │                SEC-SAXS of the complex are DEAD (dilution); no ITC
  └─ Kd > 50 µM  → Track I reduces to functional assays; redirect everything to Track D

GATE 4  N-BLOCK BINDS?  isolated N-block + gp5·TrxA
  ├─ binds (within 3-5x of full-length) → P3 confirmed; N-block becomes THE target
  └─ folded but does not bind          → P3 refuted; all interface work moves to
                                          full-length aRES; Track D proceeds anyway
```

**The important structural feature of this plan: Gate 4 can refute the project's central hypothesis without killing the rotation.** A folded, zinc-loaded, monodisperse N-block that does not bind gp5 is still the first PF18870 sample in existence and still crystallises.

---

## 5. Week 1 — exactly what to do

Nothing here needs protein, beamtime, or training. Do all of it.

**Dry (2–3 days total, not the 2–3 days per item the proposals claimed):**

1. **Cofold interface audit (P1).** Run `analysis/scripts/interface_contacts.py` across all AF3 cofold models with `--domain A:1-120:NTD1 --roi B:100-125:loop`. Report ipTM, interface PAE at the D104–S122 pair, and the per-seed contact frequency. *The cofold coordinates are not yet in `analysis/afold/` — deposit them first.*
2. **Reconcile the UC5 frame** (the I234/S293 discrepancy in §3e), then build the UC5→S_34AVR map by **two** routes: full-length MAFFT/MUSCLE MSA for the ribbon (1–48, which `hmmalign` against PF18870 **cannot** align, because PF18870 starts at 49) and `hmmalign` for 49–136, cross-checked by US-align superposition. Commit as a versioned TSV. Do **not** use an S_34AVR+gp5 cofold as a mapping vote — that is a prediction validating a prediction.
3. **Fix the Pfam accession**: the RES NADase domain is **PF08808**, not PF08808 (which is a histone-like DNA-binding family). Any HMM search built on the wrong number returns the wrong family.
4. **Specificity-determining-position scan** across T7-defending homologs (UC5, S_34AVR) vs **GN04826** (near-identical phage range, fold defence **1** against T7). Positions differing in the N-block are candidate epitope contacts. This nominates both mutants and chimera boundaries, and it is the highest-information dry item in the project — higher than the numbering map.
5. **Pre-register the crosslink prediction** (chain A = aRES, chain B = gp5) with the acidic set from §3b. Commit before the experiment.
6. **Write the cryo-EM simulation control now**, on simulated particles from the AlphaFold model at 2–3 known hinge angles with matched CTF/dose/noise: *at 143 kDa, N particles and realistic SNR, is a 13 kDa hinge rotation recoverable at all?* This is a self-contained computational project a rotation student can finish, it costs no protein and no beamtime, and a clean "no" is a genuine reportable result that saves a postdoc three months.

**Orders out (the 2–3 week synthesis clock is the real constraint):**

7. **aRES boundary series, 5 constructs, one batch** — S_34AVR 1-124, 1-136, 49-136, 40-124, plus UC5 1-136 as the one cheap bet (full-length UC5 is insoluble; isolated domains often behave where the parent did not, and the escape phages were isolated against UC5). Fold **3–4 extra PF18870 orthologs at the lead boundary into the same order** — homolog expansion is the highest-yield crystallisation rescue and ordering it now makes it available in month 2 instead of month 5. All **C-terminal TEV-His6**: C8 is residue 8, so an N-terminal tag sits on the first CxxC during folding, and six histidines compete with four cysteines for Zn(II). Cut GN04826 from the *solubility* panel and keep it as the **functional negative control** (§8.2). Drop the SUMO hedge until C-terminal has actually failed.
8. **gp5·TrxA artificial operon**, following Chiu/Tillett/March 2006 (up to 90% of gp5 soluble with TrxA co-expression vs insoluble without; 97% pure in one Ni step; ~1:1 TrxA; full activity, no exogenous Trx). **One change, and it is the opposite of what was proposed: tag gp5 C-terminally and leave TrxA untagged.** A His tag on TrxA makes Ni select for TrxA — free overexpressed 12 kDa TrxA-His6 dominates the eluate, and gp5 bound to endogenous untagged TrxA is lost in the flowthrough. The gp5 C-terminus is in the polymerase domain, far from both the TrxA-binding thumb loop (~258–333) and the epitope.
9. **gp5 variants**: S122P, D106G, **S122A** (removes the hydroxyl without proline's backbone rigidity — the only way to separate "direct contact" from "loop conformation"), and **Δ118-123** (the published 1T7P/1TK0 deletion, catalytically competent, removes G121/S122 outright).
10. **gp5 exonuclease fragment at 1–187**, with one alternative at ~1–200. **Not 1–250**: the exonuclease domain is 1–187 and the polymerisation domain is 188–698, so 1–250 ends mid-fold. Note the epitope is *inside* a T7-specific helix-loop-helix insertion (~102–122), not adjacent to one, so 1–187 keeps it intact.
11. **Three peptides**: gp5 98–132 WT, S122P, D106G. Cheap, arrive in milligrams, and they gate a lot (§8.4).

**Emails, all in week 1 — these have the longest lead times in the plan:**

12. **Sorek lab**: request the published S_34AVR aRES C-terminal-6xHis plasmid and their gp5 construct. Removes 3 weeks of synthesis, gives a known-working positive control, and resolves the methods inconsistency (pBAD is arabinose-inducible, yet 0.5 mM IPTG is specified — do not guess). Order synthetic genes in parallel so a slow reply costs nothing.
13. **Chiu/Tillett/March**: request the gp5+trxA operon plasmid.
14. **SPC** (`spc@embl-hamburg.de`): iLAB account + training. Training is mandatory before instrument access and is **on the critical path** — start it before protein exists, not after.
15. **PETRA III P64** (XAFS): ask whether Zn K-edge XANES/EXAFS on a frozen protein solution is feasible and on what route.
16. **P14** (not P13 — see §12): confirm 9.66 keV is deliverable on the tunable main endstation (the T-REXX endstation is itself fixed at 12.7 keV).
17. **EMBL Hamburg User Office**: settle the beamtime access route *this week*. The April–August 2026 external call is closed and the next call (April–December 2027) had a 26 October 2026 deadline. **A rotation starting now has no external-route beamtime inside its span** and must collect on the host group's in-house or BAG allocation.
18. **CSSB cryo-EM** (`cryoem@cssb-hamburg.de`): initial consultation, and ask in writing for (a) the current detector and energy filter on the Krios, (b) whether a rotation student counts as internal, (c) real movies-per-day.
19. **EMBL Proteomics Core Facility** (`pcf@embl.de`): confirm internal-user status, and **ask specifically whether DMTMM is available** alongside the listed BS3/DSS/DSG/EDC/DSPP.
20. **Host group**: confirm there are shakers, a lysis device and FPLC capacity for 2–4 L E. coli preps. **EMBL Hamburg SPC does not provide protein purification as a service** — it offers shared resources (an ÄKTA in a cold room; bring your own columns and loops). If the host group lacks this, the whole plan needs re-planning around PEPcore (which is at **Heidelberg**, not Hamburg) or a collaborating lab.

---

## 6. The critical path

| | | |
|---|---|---|
| **Wk 1** | All of §5 | Dry work + orders + emails. Facility training starts now. |
| **Wk 2–4** | Small-scale parallel expression screen (5 aRES boundaries × 2 temperatures ± 50 µM Zn, one strain: LOBSTR-BL21(DE3)-RIL, 16 °C overnight as published). Read on lysate SDS-PAGE + 96-well IMAC. gp5·TrxA co-expression in parallel. | **Explicitly exclude SHuffle/Origami-type oxidising-cytoplasm strains** — an engineered oxidising cytoplasm promotes disulfides among six cysteines and destroys the ribbon. |
| **Wk 4–6** | Scale the 2 best to 1–2 L. **Gate 2**: nanoDSF + far-UV CD + DLS + ICP-MS. One-afternoon limited proteolysis (trypsin only, 1:2000, 4 °C, four timepoints) as a boundary check. | Target **2–3 mg per construct for two constructs**, not 10 mg for six. |
| **Wk 5–7** | **Gate 3 — Kd.** nanoDSF FoldAffinity first (~10 µL at 0.2 mg/mL per point, no label, no surface — the cheapest Kd available and it uses the instrument you are already trained on). Then MST at two label positions. | Branch here. See §9. |
| **Wk 6–8** | **Gate 4 — does the N-block bind?** And in parallel the functional battery on the established assays (§8.2–8.4): GN04826, Δ118-123, does the peptide activate, does the N-block compete in trans, is TrxA required. | These use assays that already work and need no facility queue. |
| **Wk 7–10** | **Track D**: crystallisation plates (2 constructs × 6 screens × 2 temperatures = 24 plates, ~3.5 mg; hold the rest for optimisation). Plus one in-situ proteolysis plate — 0.15 mg for a real shot at a crystallisable core is the best-value well in the campaign. **Track I**: crosslink, quench, ship to PCF. | PCF turnaround is **4–6 weeks from receipt** — ship by week 8 or the data arrive after the rotation. |
| **Wk 9–12** | One batch-SAXS session (§8.5). Crystal optimisation / data collection if hits appear. XL-MS interpretation. | |
| **Wk 12–14** | Write-up. Hand over. | |

**Stretch goals, named as such so nobody is disappointed:** a refined apo N-block structure; a 15N HSQC; an XL-MS restraint set with full escape-mutant controls. **Not rotation goals:** any cryo-EM dataset; HDX-MS; a co-crystal; ITC; a nanobody.

---

## 7. Track D — the domain (first PF18870 structure)

### D1. Construct series and the boundary question
Five constructs (§5.7). The boundary dispute matters: the AF3 PAE block is **1–120**, the Pfam call is **49–136**, and they disagree because the N-block is bimodular — **PF18870's annotated start at 49 is exactly where the Zn ribbon ends**, so Pfam is blind to the ribbon rather than reporting a folding boundary.

Trust the PAE blocks, for one checkable reason: **R182 and R184 are defence-null and lie outside the PF08808/PF08808 call of 199–363.** A "RES domain" construct built on the Pfam boundary would physically lack two active-site arginines and be dead for the wrong reason. The PAE block 126–447 contains them.

But note the self-inflicted trap: **1–120 truncates the Pfam domain by 16 residues** and may sever it. Hence 1-124 and 1-136 in the series, and hence limited proteolysis as an empirical arbiter. Also fix any full-length construct at **1–447**, not 1–437 — the PAE block ends at 447 and deleting 438–447 from an ordered block has no justification.

**Means:** 49-136 soluble while all ribbon-containing boundaries fail → the ribbon needs RES-side context, the N-block is not autonomously bimodular, and Zn-SAD/XAS is unavailable. 1-120 soluble while 49-120 is not → direct evidence the ribbon is required for NTD-1 folding, which vindicates the PAE boundary over Pfam and is a real result about the family.

### D2. Fold and metal — the gate that must not be soft
- **nanoDSF** (Prometheus NT.48, 48 capillaries, ~10 µL at 0.2–0.5 mg/mL): buffer screen as **two runs of 48**, not one of 96. Cooperative transition, Tm > 45 °C.
- **Far-UV CD is mandatory, not optional.** A cooperative Tm is a weak fold gate for a family with no experimental structure anywhere — a compact non-native species gives one too. NTD-1 is predicted mostly α-helical, so a CD spectrum is a direct test of the prediction and the cheapest publishable statement anyone has made about PF18870. (10 mM HEPES pH 7.5 + 50 mM NaF; no Tris, no phosphate, TCEP ≤1 mM.)
- **Zinc, by measurement rather than inference.** Primary: **ICP-MS or ICP-OES** (~0.2 mg, absolute, definitive). Secondary: **native ESI-MS** — holo vs apo differ by **~+63.4 Da**, not +65.4, because Zn²⁺ displaces two thiol protons; it needs micrograms, works perfectly at 14 kDa (where mass photometry does not), and reports oligomeric state at the same time. Tertiary: PAR, *but* PAR is a general divalent-metal indicator (Co, Ni, Cu, Cd, Hg, Mn) run on Ni-NTA-purified protein, so it needs a protein-free blank, a Ni²⁺ standard curve alongside the Zn²⁺ one, and a reductant-free aliquot. **Do not run PMPS on protein containing TCEP** — organomercurial zinc release is reversed by reductants.
- **Chelator challenge, interpreted asymmetrically.** Use **TPEN** (apparent Zn Kd ~6.4×10⁻¹⁶ M), not EDTA (~2.3×10⁻¹⁴ M). A large Tm drop is informative; **a flat response is uninformative and must not be read as absence of metal** — a folded Zn(Cys)₄ site (Kd 10⁻¹³–10⁻¹⁷ M, slow exchange) is thermodynamically and kinetically resistant to EDTA. One proposal planned to abandon the whole zinc strategy on exactly this null artefact.
- **Zn K-edge XANES/EXAFS at PETRA III P64** if access allows. On the same campus, directly measures coordination number, ligand identity (S vs N/O) and Zn–S distance on a frozen solution sample — strictly more than Co(II) UV-vis can report, on a site nobody has ever observed. This is the best single idea to come out of the critiques.
- **Co(II) substitution**: only if you budget properly. d–d bands are ~10²–10³ M⁻¹cm⁻¹, so 200 µg of a 14 kDa protein in 1 mL gives A ≈ 0.005–0.015 — noise. Needs ~0.5 mg in a 100 µL microcuvette. Otherwise skip it in favour of XAS.
- **Buffer discipline, correctly weighted.** The real hazards are **failed metalation during cytoplasmic folding** and **cysteine oxidation** (six cysteines, §3d) — not EDTA stripping a folded site. So: 50–100 µM ZnSO₄ in the medium at induction; TCEP at 1–2 mM; degassed buffers; avoid phosphate (Zn precipitates). Avoid EDTA as cheap insurance, but do not build the prep around avoiding imidazole — 250–500 mM imidazole cannot strip a genuine structural tetrathiolate, and on-resin TEV cleavage (itself unvalidated in 25 µM Zn / low TCEP) adds a failure mode to solve a non-problem. **Desalt into zinc-free buffer before any zinc quantification** — 25 µM free Zn against 35–140 µM protein is a 20–70% background that makes the stoichiometry meaningless.

### D3. Crystallisation
24 plates, two constructs, 2 concentrations (10 and 20 mg/mL), two temperatures (RockImager 1000 at 19 °C, RockImager 182 at 4 °C), monitored remotely via CRIMS. Plus one in-situ proteolysis plate.

**Delete every reference to CrystalDirect.** There is no CrystalDirect plate system or harvester at EMBL Hamburg — it was developed at the Grenoble Outstation and the harvester sits on MASSIF-1 at the ESRF. Plan on manual harvesting from standard SWISSCI/MRC sitting-drop plates, and confirm the stocked screen catalogue by email rather than assuming six named screens.

Be realistic: a 14 kDa all-α domain tethered to a zinc ribbon across a flexible hinge is **not** a favourable crystallisation shape — high solvent content, few flat surfaces. ~35–45% chance of any hit.

**Rescue, if dry:** homolog expansion only (already in the week-1 gene order, so available from month 2), with the original construct always on a side-by-side control plate so a hit is attributable. Surface-entropy reduction is second, and its exclusion zone must be **the whole candidate epitope face** defined from the UC5/GN04826 divergence pattern, not just a few residues around D104. Reductive methylation and microseeding are afterthoughts. **Nanobodies/sybodies are out** — 2–3 months before a candidate exists, and you then still have to crystallise the complex.

### D4. Phasing — inverted from what was proposed
**Collect the best-resolution native dataset you can and solve by molecular replacement with the AlphaFold model.** At mean pLDDT ~90 Phaser will probably solve it in minutes. "No close homolog, so MR is unavailable" is no longer true.

Then add a **short low-dose sweep at the Zn peak purely for the anomalous difference Fourier** — ~30 extra minutes, and it still delivers experimental proof of the Zn site plus an independent check on sequence register.

Do **not** design the collection around Zn-SAD. The proposed protocol (reduced transmission, anomalous multiplicity 20–40, merging across 3–10 crystals) trades resolution for redundancy, and for a structure MR will solve outright, **resolution is the deliverable**. The risk is spending a single precious shift on a 2.8 Å high-multiplicity dataset when the same crystals would have given 2.2 Å.

Two corrections worth carrying:
- **The Bijvoet ratio is 2.69%, not 1.90%.** The canonical Hendrickson/Smith estimate is `sqrt(N_A/2N_T) × (2 f″/Z_eff)`; one proposal used `√2 f″/Z_eff` where the published form has `2 f″/Z_eff`. I recomputed: 1 Zn, 120 residues (N_T ≈ 936), Z_eff 6.7, f″ 3.9 e → **2.69%** (1-124: 2.65%; 1-136: 2.53%). The error is conservative, so the anomalous experiment is *easier* than claimed — which argues for ~10-fold multiplicity, not 20–40.
- **The 49–136 construct contains no zinc** (all four Cys are at 1–48), so any Bijvoet ratio quoted for it is meaningless. It is an MR-only target.
- **SeMet is a live hedge** (§3c): five Met in 120 residues. Run one SeMet trial expression if zinc proves substoichiometric.
- An apo N-block structure **cannot** report on hinge mobility — every ribbon-containing construct terminates at 120, 124 or 136 and cuts the hinge (120–125) out. Do not promise "quantitative B-factors for the hinge" as a deliverable.

### D5. The domain's solution deliverable, which needs no crystal
**Batch-mode SAXS of apo full-length aRES against the AlphaFold model.** 50 µL at 1–3 mg/mL (~150 µg), three concentrations, no partner, no Kd, no complex, no gp5. It is the only structural experiment in the whole plan immune to the material and affinity gates, and it tests the architecture everything else assumes: the repo model gives Rg 31.9 Å with the two block centroids 50.5 Å apart, where a compact 52 kDa globule would be ~25 Å. A CRYSOL fit plus the Kratky shape either confirms that extended two-block arrangement in solution or shows the modelled hinge geometry is wrong.

**Also measure the sensitivity before committing:** on aRES alone the N-block is 27% of the scattering mass and a 5 Å change in block separation shifts Rg by ~1.6 Å (~5%, measurable). In the 143 kDa complex the N-block is ~10% of the mass and the same change shifts Rg by ~0.6 Å (~1.7%, at or below the floor). **SAXS is well-posed on apo aRES and ill-posed on the complex** — which is the reverse of how it was staged.

---

## 8. Track I — the interaction

### I1. Kd first. Everything else is downstream of one number.
No affinity number exists anywhere in the literature. The paper's evidence for direct binding is an in-cell Ni-NTA pull-down from lysate plus an in vitro **activity** reconstitution — and an activity readout is satisfied by transient, sub-stoichiometric encounter. A sensor switched on by a brief collision and staying on gives clean ADPR-1P signal at Kd = 10 µM. "The paper showed binding" does not mean "the complex is SEC-purifiable".

Ladder, cheapest first: **nanoDSF FoldAffinity** (isothermal analysis of ligand-dependent Tm shift; ~10 µL at 0.2 mg/mL per point; monitor the aRES transition and run gp5·TrxA alone as a reference trace) → **MST** (RED-tris-NTA via the His tag, not NHS chemistry; ~1 µg labelled protein and ~0.2 mg gp5·TrxA per titration) → **ITC** only if Kd lands in 0.15–3 µM.

Two MST caveats to design around: the published C-terminal His6 sits at the end of a **disordered tail (448–459, pLDDT 54)**, so the dye may be 40+ Å from the N-block epitope and simply not sense a distal binding event — **run both a C-terminal and an N-terminal His construct and treat a null as negative only when both are flat**. And keep TCEP at 0.25 mM for MST plates; it reversibly quenches Cy5-class dyes.

### I2. The controls that make everything else mean something
These appeared in the "missing experiments" list of **all four** expert areas. They are cheap, they use assays that already work, and the plan has no false-positive detection without them.

- **GN04826 aRES as the natural non-binder.** Near-identical phage range to UC5 but **fold defence 1 against T7**. A correctly folded, independently evolved protein that should not respond to the trigger — immune to the "your mutant is misfolded" objection in a way no point mutant is. Run it alongside WT in every binding, crosslinking, footprinting and activation experiment.
- **A non-trigger ligand**: S_34AVR defends T7, SECφ18 and SECφ4 but not T5 or T4, so T5 DNA polymerase or T4 gp43 should give no binding. Without it, every Kd is uncontrolled for the sticky-protein artefacts that are the characteristic failure mode of SPR and MST at micromolar analyte. TrxA alone is far too easy a test.
- **gp5 escape and alanine variants with a folding gate.** S122P and D106G are the phenotype; **S122A** separates contact from conformation; **Δ118-123** is the published, catalytically competent loop deletion. Note the reasoning on D106 must be the right way round: **D106G was isolated as a viable escape phage, so that polymerase folds and replicates a genome.** The 0.6% RSA burial argument predicts misfolding for the *novel* D106A/D106N, not for the escape mutation itself. Gate the novel ones by Tm + primer extension on the reconstituted 1:1 complex with 5 µM free TrxA.
- **aRES D104A/D104N — removed from the week-1 gene order.** Building point mutants on a single
  low-ipTM cofold pose is the clearest case of over-committing to an unvalidated prediction. Order
  them only once the docking runs or the GN04826 divergence analysis nominate a surface
  independently. The phenotype has also never been shown.
- **A one-plate "is anything missing?" control**: activation ± DNA, ± dNTPs, ± Mg²⁺, ± aRES-null E. coli lysate. One day, and it protects every negative result in the project from the "you were missing a cofactor" objection.

### I3. The interface restraint set — XL-MS, with the chemistry fixed
This is Track I's realistic deliverable. ~100 µg, no crystals, no beamtime, no microscope queue, 4–6 weeks at the EMBL PCF (BS3/DSS/DSG/EDC/DSPP listed; 50 µg per sample at ≥0.2 mg/mL; user does the chemistry to the quench step).

**Chemistry, decided by §3a:**
- **EDC, or preferably DMTMM, is the sole interface chemistry.** Zero-length Lys↔Asp/Glu matches this interface because the predicted contact is to an acidic residue and there is a real acidic patch around D104. DMTMM works at pH 7.5; EDC needs pH ~6 in MES, activates every carboxylate in a 143 kDa complex, is known for non-specific laddering at ~7 µM protein, and asks the complex to survive pH 6 — untested. Ask PCF for DMTMM.
- **Run it as discovery, not confirmation.** Search the whole acidic face of the N-block, not a
  pre-chosen pair. *Write the prediction down first* — gp5 K103/K114/K118 → aRES D104, D106, D97,
  E98 or D112 — so a hit on it means something; but do not design the experiment around it, because
  it comes from one low-confidence pose. A crosslink anywhere in 49–120 answers P3.
- **BS3/DSS is demoted to the intramolecular folding control** that validates the S122P and D106G preps. It cannot report NTD-1 and a RES-block hit does not refute anything (§3a).
- Sulfo-SDA / genetically encoded pBpa would capture a Ser–Asp contact directly but need bespoke non-specific database searching and are not rotation-scale.
- **A crosslinked sample can never feed a conformational-heterogeneity analysis.** Always keep an uncrosslinked reference.

**Means:** crosslinks placing the gp5 100–125 loop against aRES 49–120, lost in S122P/D106G while the intramolecular maps are unchanged, present in WT and absent in GN04826 → the escape genetics are physically linked to a defined interface. That is the strongest result available in this area without a structure, and the restraints then fit directly into any low-resolution map. Crosslinks to the Zn ribbon instead → the ribbon, not NTD-1, is the contact surface, which is a different and interesting answer.

### I4. Functional tests that cost almost nothing and gate a lot
All on the established NAD/NADH-Glo and LC-MS ADPR-1P assays. **Run these in weeks 6–8, not at the end.**

1. **Does the gp5 98–132 peptide activate aRES?** Direct functional test of whether the epitope loop is the trigger. One peptide, one plate. Gates all peptide-based NMR and the sensitive version of the SAXS experiment.
2. **Does the isolated N-block compete in trans?** Excess N-block should inhibit gp5-driven activation of full-length aRES. Uses only constructs already required, on an assay that already works — the strongest single functional test of P3.
3. **Is TrxA required?** gp5 alone (needs a trxA-deficient host) + aRES + NAD⁺. Either removes TrxA from every downstream buffer and construct decision, or makes it part of the recognised surface. If gp5·TrxA is required for recognition, aRES is reading a **host-factor-captured replication machine** — which would make the λ puzzle (λ encodes no DNA polymerase, yet two aRES proteins defend against it; its only replication proteins are O and P) much more tractable.
4. **Does gp5(Δ118-123)·TrxA still activate?** A clean loop-deletion control on an already-validated construct, and a far stronger test of the epitope than any map at 6–9 Å.
5. **Does aRES inhibit polymerase activity?** Primer extension or dNTP incorporation. Tests a second possible abortive-infection mechanism in one day, with a readout nobody can dismiss as a gas-phase artefact.

### I5. Protease protection — keep the good half
**Keep the aRES arm.** aRES is a two-block protein with a flexible hinge at 120–125, so a limited trypsin or elastase digest should cut there and give two gel-resolvable fragments of ~14 and ~38 kDa. That makes it simultaneously (a) a test of whether gp5 binding opens or closes the hinge — *new* cleavage on binding is the autoinhibition-release signature — and (b) the classical limited-proteolysis **boundary screen**, which is worth more to Track D than the protection mapping itself. It also works at weak affinity and needs no instrument booking.

**Drop the gp5-loop protection arm.** gp5 is 704 residues with ~80 trypsin sites; losing three (K103/K114/K118) does not produce a visible change on a 4–20% gel, and aRES and gp5 fragments overlap across the whole 10–60 kDa range with no discriminator in the design.

### I6. Sample QC, and where it belongs
Do it at the SPC, not at a beamline. Analytical SEC + DLS + nanoDSF + A280 answers "monodisperse, right mass, survives a column" on 20 µg in an afternoon with no proposal cycle. Spending a half-day of P12 on quality control is a category error.

**Add two gates the QC lists omit:**
- **A260/A280 on every gp5 prep** (and on aRES). gp5 is a DNA polymerase purified from E. coli; co-purifying nucleic acid inflates A280-derived concentration, dominates small-angle scattering, and corrupts MALLS mass and mass-photometry masses with no obvious warning. Include a heparin or Q step and nuclease treatment.
- **NAD/NADH-Glo on every WT aRES batch** before it goes into a capillary or a plate — and confirm samples are NAD⁺-free before any hours-long incubation of WT aRES with its activator.

**Mass photometry: know its limits.** The SPC instrument is listed as a **Refeyn One** on one source, not a TwoMP, so the "30 kDa–5 MDa" figures may not apply. Either way: it **cannot measure the 14 kDa N-block** (below every quoted floor of 30–50 kDa), 52 kDa aRES sits at or below the vendor-calibrated range, and with ~25 kDa FWHM at 66 kDa the 52 and 91 kDa peaks overlap substantially. Use it as **prep-quality QC on full-length aRES and gp5·TrxA only**. It also works at 5–20 nM, so **absence of a 143 kDa peak is not evidence against the interaction** — any complex with Kd ≥1 µM is >98% dissociated there. The one genuinely good use: measure the **quenched, crosslinked** sample, where the covalent complex is immune to dilution and you read out crosslinking efficiency before paying the MS invoice.

---

## 9. Decision points

| # | Trigger | If yes | If no |
|---|---|---|---|
| **1** | Does **any** NTD-1 contact reproduce across cofold seeds and across the two docking runs, and at what ipTM? *(Not: does D104–S122 specifically reproduce)* | The N-block is the contacting domain on two independent methods → interface work can be aimed | No convergent NTD-1 contact → do **not** order interface mutants. Track I continues as Kd + isolated-N-block binding + trans competition, all of which are residue-agnostic. Track D is unaffected |
| **2** | aRES monomeric by SEC-MALS/AUC? | 1:1 stoichiometry arithmetic holds | Re-do every mass, mask and particle-count figure. A ~100 kDa aRES dimer is unresolvable from 92 kDa gp5·TrxA by mass photometry (6GW6 precedent) |
| **3** | N-block soluble + cooperative Tm + CD-consistent + ~1 Zn? | **Track D opens.** Crystallisation plates, SAXS, (NMR) | Boundary redesign from limited proteolysis; if all boundaries fail, the N-block is not autonomously folding — a real result about the family — and all structural work moves to full-length aRES |
| **4** | gp5·TrxA ≥0.5 mg by week 6? | Track I opens | Track I reduces to functional assays on the commercial enzyme. Put everything into Track D. **Decide this honestly at week 6, not week 10** |
| **5** | Kd < 1 µM? | Everything works, including complex SEC-SAXS and low-concentration crosslinking | 1–20 µM: keep XL-MS + protease protection + functional assays; kill mass photometry of the complex, complex SEC-SAXS and ITC. >50 µM: Track I is functional-only, and suspect the pull-down reported avidity or a lysate cofactor |
| **6** | Isolated N-block binds gp5·TrxA within 3–5× of full-length? | **P3 confirmed.** The N-block becomes the structural target for everything, and a 14 kDa domain + small partner is a realistic co-crystal | **P3 refuted** (provided Gate 3 passed — this is why the fold gate must precede the binding test). The interface needs the RES block or the hinge; the cofold is wrong about which block contacts gp5; all interface work uses full-length aRES. Track D is unaffected |
| **7** | Crystals by week 10? | Optimise, collect, solve by MR + Zn anomalous validation | Homolog expansion (already in hand from the week-1 order); and the apo-aRES SAXS + CD + ICP-MS + proteolysis-defined boundary package is still a coherent first characterisation of PF18870 |
| **8** | Does the peptide activate / does the N-block compete in trans? | Both positive → a small, crystallisable, functionally validated minimal system. This is the best possible outcome of the rotation | Peptide-negative is weakly interpretable (§11); competition-negative is strong evidence against the isolated N-block being the sensor |

---

## 10. Dropped, and why — one line each

So you know these were considered and rejected rather than overlooked.

**Cryo-EM (the whole arm, as a rotation deliverable).** The informative N-block is 13 kDa — 9% of a 143 kDa particle — on a pLDDT-54 hinge; the most likely successful outcome is a 3.5 Å map of gp5·TrxA + RES domain with the answer-bearing lobe smeared into a 10 Å cloud. Protein by week 4 + screening + collection + processing lands at week 14–15 with zero slack. **Keep**: the gp5·TrxA-alone same-session control and the post-dock unassigned-volume check (the two defences against silent dissociation producing a beautiful map of nothing), the decision record that gp5·TrxA already *is* a 92 kDa rigid fiducial with a 2.2 Å reference so no nanobody/megabody/MBP fusion is needed, and the simulation control run on simulated data in week 2. **Note**: cryo-EM cannot resolve the S122–D104 hydrogen bond — that loop is deleted from 1T7P and exists only in models. If that contact is the actual goal, cryo-EM is the wrong instrument.

**Focused classification of the N-block.** Gated on a Krios dataset that will not exist; its negative is interpretable only if a multi-week matched-SNR simulation control passes; and the apo reference it needs for the autoinhibition comparison does not exist.

**ITC.** For c = 10 at a 20 µM Kd you need ~200 µM gp5·TrxA in the cell — **~5.5 mg of a 92 kDa complex per titration**, an order of magnitude beyond supply, in the weak-binding scenario the plan itself calls plausible. Cut it in week 1 so that protein goes into the activation dose-response instead. (Its protonation-change extension was also justified on "a His/Asp-containing interface" — no histidine has been placed at this interface anywhere.)

**HDX-MS.** Blind to the hinge *by construction* — its own internal controls specify that 120–125 must exchange to completion in the first timepoint, and a region fully exchanged at 10 s cannot report a bound-vs-apo difference. Plus 2–4 mg of protein and no verified local access. If pursued at all, aim it at a route that demonstrably exists (the Astbury facility at Leeds is an Instruct-ERIC access platform with a documented HDX pipeline) and treat it as handover.

**gp5 exonuclease fragment as a co-crystallisation partner.** In Pol I family polymerases the proofreading exonuclease domain is structurally intertwined with the polymerase domain rather than tethered; no literature precedent for it folding alone; and it is the load-bearing assumption of the whole co-crystal rescue path.

**The full aRES·gp5·TrxA co-crystal.** ~145–160 kDa with a flexible hinge, a flexible 392–408 loop, a disordered C-tail and an epitope loop never seen experimentally. Will not crystallise in a rotation. Do not plan beamtime around it.

**Native MS of the ternary complex for TrxA retention.** gp5·TrxA has Kd ~5 nM, so retention is the expected result and tells you nothing, while the interesting outcome (aRES ejects the processivity factor) is exactly what cannot be distinguished from desolvation-driven dissociation. Replaced by the polymerase-inhibition assay (§8.4.5). *Native MS for the zinc question is kept* (§7.2).

**The Cys4→Ala/Ser quadruple mutant.** Almost certainly unfolded, so "no zinc and no fold" cannot distinguish "zinc is required" from "this mutant is broken". Keep single C→A/C→S variants, where a retained fold with reduced metal affinity *is* interpretable, and budget for them failing.

**GraFix.** Its glycerol must be removed before freezing, and that removal is itself reported to damage fragile complexes.

**Synthetic peptide as a *binding* reagent.** A 26–35-mer excised from a helix-loop-helix will be random coil, and the working model already says D106 is a buried intramolecular scaffold (0.6% RSA) holding the loop — a free peptide has no scaffold. **But buy the peptides anyway** and use them *functionally* (does the peptide activate?) and as an ITC/MST probe against full-length aRES: the positive is strong and cheap, and it hands crystallography a 14 kDa + 35-mer target. Just do that test *before* committing three months of NMR negotiation and 1 mg of labelled protein.

**NMR Tier 2 (peptide CSP titration).** Gated on the above. Tier 1 — one HSQC per boundary construct as a fold/Zn check, with 50–100 µM ZnSO₄ in the 15N M9 (minimal medium is zinc-poor and an unloaded ribbon will not fold) — stays as a stretch goal. Fix the dynamics design to per-residue 15N R2 and {1H}-15N hetNOE compared between 1-120 and 49-136; "two tethered modules give an anomalously short τc" is not correct, a tethered pair tumbles with an *intermediate* τc.

**P13 for anything zinc.** Its energy is listed as **fixed at 12.7 keV** as of November 2025 due to technical issues; the Zn K-edge is 9.659 keV. Use P14. One consolation: XRF detection of Zn Kα only needs incident energy *above* 9.659 keV, so the metalation check survives even at 12.7 keV.

**CrystalDirect.** Not at Hamburg (§7.3).

**Cell-free translation as a toxicity bypass.** Cell-free lysates contain NAD⁺, so an active NADase depletes the reaction.

**pBAD-vs-IPTG guessing.** The published methods are internally inconsistent (pBAD is arabinose-inducible, yet 0.5 mM IPTG is specified). Request the plasmid; do not build on a guess.

**Over-scoped versions of things kept:** the 12-homolog × full-boundary matrix (350–450 cultures, not "24 plasmids / 6–8 days") cut to 5 constructs × 3 conditions; the 25–30-condition cryo grid screen cut to a 6–8 grid Round 1; the "10 mg per construct × 6 constructs" budget cut to 2–3 mg × 2; the 96-condition nanoDSF screen split into two 48-capillary runs; the three-protease six-timepoint proteolysis matrix cut to one afternoon of trypsin.

---

## 11. Experiments whose negative result means nothing — flagged

Run these knowing the asymmetry, and never report the null as a finding.

| Experiment | Why the negative is uninterpretable |
|---|---|
| **Mass photometry of the complex** | Measures at 5–20 nM; any complex with Kd ≥1 µM is >98% dissociated. Absence of a 143 kDa peak says only "not sub-nanomolar" |
| **Negative-stain EM** | Uranyl acetate and carbon routinely destroy complexes that are fine in vitreous ice. A field of aggregates does not license "the sample is bad" — only when solution-state *and* stain both look bad |
| **BS3/DSS interface crosslinking** | NTD-1 has no lysine (§3a). The chemistry is structurally incapable of the readout, so its null carries no information. And a RES-block hit at 29–32 Å does not refute NTD-1 |
| **EDTA thermal shift as a metalation test** | A folded Zn(Cys)₄ site resists EDTA thermodynamically and kinetically. Flat ≠ apo. Use TPEN, and settle it by ICP-MS |
| **Flat SPR sensorgram** | Immobilising aRES may occlude the N-block epitope. Confirm in free solution before reporting |
| **Flat MST trace** | The dye may sit on the disordered 448–459 tail, 40+ Å from the epitope. Requires both label positions flat |
| **Peptide CSP / peptide binding negative** | An isolated loop out of context may simply not adopt the presented conformation |
| **gp5-loop protease protection** | Three sites out of ~80; gel-invisible. (Dropped for this reason) |
| **Failure to crystallise** | Says nothing about whether the domain is folded or binds gp5. Interpretable only alongside the CD/nanoDSF fold data |
| **Cryo-EM failure to see the N-block** | Says nothing about whether the N-block is the sensor in solution. Must not be reported as evidence against the cofold |
| **Cryo-EM hinge-motion null** | Meaningless unless the matched-SNR simulation control proved the pipeline *could* have recovered a known rotation |
| **SAXS: no Rg change on binding** | Excludes only a *large* change in radial mass distribution, and in the complex the N-block is ~10% of the scattering mass. Report as "no large compaction or extension detected", nothing stronger |
| **Mutant binding loss without a folding gate** | Indistinguishable from misfolding. Always pair with Tm + unchanged intramolecular crosslink map + (for gp5) primer extension on the 1:1 complex |

---

## 12. Facility claims to verify on arrival — the single list

None of the five experts could read `embl.org`, `embl-hamburg.de` or `cssb-hamburg.de` directly (403 through the session proxy), so **everything below is from search snippets and must be confirmed locally.** Several are load-bearing.

**Settle in week 1 (critical path):**
1. **Beamtime access route** for a rotation student — external calls are closed; confirm the host group's in-house/BAG allocation with the User Office.
2. **Does the host group have shakers, a lysis device and FPLC for 2–4 L preps?** SPC does *not* offer protein purification as a service — only shared resources (ÄKTA in cold room 48e, bring your own columns/loops; Beckman Optima XPN-80). No lysis device appears on any list.
3. **SPC training and iLAB queue lead time**, and the fee structure for an internal rotation student.
4. **P14**: is 9.66 keV deliverable on the tunable main endstation? (listed 7–20/35 keV, EIGER2 CdTe 16M, AMPTEK SDD for fluorescence; T-REXX is fixed at 12.7 keV). And is a fluorescence energy scan a documented, supported user workflow? Contacts named in the material: Gleb Bourenkov, Thomas R. Schneider.
5. **P13's actual state** — reported fixed at 12.7 keV as of Nov 2025; also sources disagree on its detector (PILATUS 6M-F vs EIGER 16M) and upper energy limit (17.5 vs 18 keV).

**Instruments — confirm presence and model:**
6. **Which Refeyn?** One source says Refeyn **One** (MOSBRI EMBL-SPC page), not TwoMP. The models differ materially in mass range and working concentration, and the whole stoichiometry argument depends on it.
7. **SEC-MALS**: does not appear on any SPC instrument list anyone could retrieve. The only verified MALLS is the in-line detector train on P12 (MALLS + UV + RI + QELS). This matters because mass photometry cannot cover the small constructs — so if SEC-MALS is absent, native MS or analytical SEC against standards is the substitute for stoichiometry, and the aRES oligomeric-state question (P5) needs a plan.
8. **CD spectrometer at Hamburg?** The "molecular-biophysics" page listing CD, MALDI with high-mass detector, FTIR, SPR and Octet is titled `molecular-biophysics-hamburg-embl` in one retrieval but was attributed to the Heidelberg PEP core in another. **Far-UV CD is a mandatory fold gate in this plan**, so resolve it. One email.
9. **Biacore T200 vs Octet RED96** — both appear on CSSB/EMBL-SPC listings, but the T200 may sit in the Chemical Biology core. Prefer the **Octet** anyway: far less method development than a Biacore surface for a rotation timescale. Note **KinGenie** (the eSPC BLI-kinetics module) is explicitly beta — plan to fit in the Octet software too.
10. **ICP-MS or ICP-OES** anywhere on the DESY / University of Hamburg campus. No expert found evidence either way, and this is now the *primary* zinc method. If genuinely absent, it becomes an external collaboration to arrange in week 1.
11. **Native ESI-MS** access. SPC lists MALDI-TOF, which is unsuitable for intact non-covalent metal stoichiometry. The Uetrecht group (native/top-down MS) works at CSSB on the DESY campus, with the professorship now at Lübeck — terms for a rotation student unverified.
12. **Intact-protein LC-MS** — is it user-bookable at SPC, or does band identification go to the UKE proteomics core? Affects the proteolysis boundary work.
13. **A 120 kV TEM** for negative-stain imaging. The Talos L120C appears only in a March 2026 University of Hamburg maintenance tender; CSSB's own page lists the Krios G3 and Talos Arctica. SPC prepares negative-stain grids but no TEM model at EMBL Hamburg was found.
14. **CSSB microscope inventory generally** — the K3 detector and BioQuantum filter that collection parameters were tuned to come from 2019–2020 mailing-list job adverts, not a current page. No source confirms a Krios G4, Selectris, Falcon 4i or cold FEG. Whether any Volta phase plate is operational and offered. Whether a Chameleon or any non-blotting plunge-freezer exists (plan assumes Vitrobot Mark IV only).
15. **PETRA III P64** — XAFS beamline, ~4–44 keV, 100-pixel Ge fluorescence detector, specified for dilute metalloproteins. Confirm access route and frozen-solution sample requirements.
16. **EMBL PCF (Heidelberg)** — confirm EMBL Hamburg counts as internal; confirm the quench reagent in their protocol, the search engine, fees, and **whether DMTMM is available**. Verified: BS3/DSS/DSG/EDC/DSPP, 50 µg per sample at ≥0.2 mg/mL, user does chemistry to quench, **4–6 weeks turnaround from receipt**.
17. **PEPcore construct-screening access** for a Hamburg-based rotation student, cost and turnaround — it is at **Heidelberg**, and commercial vectors there are shared with internal EMBL users only.
18. **P12 sample specs** (from snippets, to re-check): batch SAXS 10 µL minimum / 50 µL recommended / 100 µL for a concentration series, ~90 s per sample; SEC-SAXS 20–100 µL at 5–10 mg/mL, ≥500 mL mobile phase, 15 min on a 3 mL column / ~50 min on 24 mL; mail-in supported; **SEC-SAXS use must be justified in the proposal** — pure QC is not a justification.
19. **SPC crystallisation specifics**: exact screen catalogue, drop volumes, plate types (SWISSCI/MRC — *not* CrystalDirect). Verified present: Mosquito-LCP, RockImager 1000 (19 °C), RockImager 182 (4 °C), Scorpion, Liquidator96, plate sealer, CRIMS.
20. **Luminescence plate reader**: the Tecan **Spark 20M** on the SPC/CSSB list has a dedicated luminescence PMT (glow/flash/multicolour), so NAD/NADH-Glo should be fine — confirm. And campus LC-MS access for ADPR-1P.
21. **Protein NMR**: no biomolecular spectrometer found at EMBL Hamburg or CSSB. The University of Hamburg platform (300–600 MHz Bruker, coordinator Dr Thomas Hackl) is described as solution small-molecule work. Realistic route is a partner site via Instruct-ERIC (e.g. BNMRZ at TU Munich, which requires the labelling pattern agreed with their staff in advance). **Confirm a specific spectrometer and access route before labelling any protein.**
22. **PMPS / organomercurial procurement and disposal** under current mercury rules — a real question for a 12-week rotation, and a further reason to make ICP-MS primary.

**Also unverified, non-facility:** NEB M0274S's mass concentration (verified: 300 units at 10,000 U/mL = **30 µL of stock**; the "250 units" in the background is wrong; **no mg/mL is published**, and Chiu's thesis reports the commercial prep carries a **Phe102→Ser substitution in thioredoxin** and is less processive at 37 °C, so it is not even the same molecule — buy exactly one tube as an activity reference and produce everything else); gp5 and gp5·TrxA mg/L yields from a modern His-tag co-expression (the ~700 mg/12 L patent-era figures used a separate thioredoxin booster culture and are not transferable); the S_34AVR sequence and accession (all five aRES IMG scaffold IDs in `docs/04_data_inventory.md` were transcribed off the small type in Fig. 1A and are flagged VERIFY — **get the sequences from Table S4 or from the Sorek lab, not from a figure**); the canonical gp5 exonuclease-deficient substitutions (look them up, do not guess); TEV cleavage efficiency in 25 µM Zn and low TCEP (test on 100 µg first); and Co(II)-Cys4 band positions and extinction coefficients.

---

## 13. Honest scope

**What a 12–14 week rotation will plausibly produce:** the numbering map and the UC5 frame reconciliation; the cofold interface audit; the aRES oligomeric state; a 5-construct expression screen with a defined lead boundary; a fold and zinc characterisation of the N-block (nanoDSF + CD + DLS + ICP-MS, possibly XAS) — **which is already the first experimental characterisation of any PF18870 domain**; a measured Kd or a defensible lower bound; the functional control battery (GN04826, Δ118-123, peptide activation, trans competition, TrxA requirement); one batch-SAXS session on apo aRES against the model; crystallisation plates set and imaging under way; and an XL-MS sample shipped.

**What will probably be half-done:** crystal optimisation; the XL-MS interpretation (the facility clock, not the student, decides); a 15N HSQC.

**What will not happen, and should be promised to nobody:** a solved co-complex structure; any cryo-EM dataset; HDX-MS; ITC; a nanobody.

**Two scheduling facts that dominate everything.** Gene synthesis is 2–3 weeks, so orders go out in week 1 on best-guess boundaries that the alignment refines by PCR afterwards. Facility training and beamtime access are longer lead times than any experiment, so both start in week 1 — before protein exists.

**And one framing point worth saying to the PI at the start.** Crystallography of a never-before-crystallised domain family *is* mostly construct engineering; by time, the crystallography arm is roughly 60% cloning and protein production, which is the thing the student said they did not want. The honest response is not to hide that but to point at the techniques that need **no new constructs** on the one prep that already exists — SAXS at P12, SPC biophysics, and XL-MS. Those are the fastest route to a structural-biology result inside a rotation, and the apo N-block structure is the one with the highest ceiling if the protein cooperates.