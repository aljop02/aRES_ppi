# Hypothesis: PF18870 (NTD-1) is a modular phage-protein sensor domain

> **Update (Oct 12, 2026).** This document was written before the existing work in `notes/vault/`
> was available. It is still correct in its logic, but two things are now better established and
> one framing here is superseded — see **`05_status_after_vault_integration.md`**:
> - The existing aRES×gp5 cofold placed a specific residue pair across the interface. **That reading
>   has been dropped** — it was one seed at poor ipTM, it was never reproduced, and nothing in the
>   project rests on it any more ([`09` §4](09_blind_zn_run_results.md)). H1 is supported by domain
>   architecture and the escape-mutant genetics, *not* by any cofold contact.
> - **PDB 1T7P lacks residues 118–123**, so the epitope loop is not in any experimental structure.
>   The "well-defined helix-loop-helix" framing in §3 below overstates what is known.
> - Real domain boundaries: N-block **1–120**, hinge 120–125, RES side **126–447**
>   (InterPro: PF18870 49–136, PF08808 199–363).

Working document for the wet-lab rotation. Status labels used throughout:
**[FACT]** = shown in Osterman et al. 2026 or another cited paper ·
**[INFERENCE]** = my reasoning from facts, testable ·
**[SPECULATION]** = plausible, unsupported, listed to be falsified.

---

## 1. Hypothesis statement

> **H1 (narrow).** The N-terminal PF18870 domain ("NTD-1") of aRES is the module that binds the
> T7 DNA polymerase (gp5), and this binding is what activates the C-terminal RES NADase.
>
> **H2 (broad).** NTD-1 is a *modular, transplantable phage-protein sensor*. It has been
> recombined onto several unrelated toxic effectors (RES NADase, HEPN RNase, ABC ATPase), and the
> phage-recognition specificity of each system is determined by its NTD-1 sequence, not by its
> effector.

H2 is the ambitious version and the one that would make the rotation interesting beyond a single
system. H1 is the prerequisite.

### Why this is worth doing

The authors themselves flag it as unresolved **[FACT]**:

> "Presumably, the N-terminal domain of the aRES protein (Pfam accession PF18870) […] is
> responsible for recognition of the T7 DNA polymerase. […] Future studies will be necessary to
> understand the structural and molecular basis for DNA polymerase recognition by aRES and how
> this recognition activates NAD⁺ cleavage."

So: not a crowded field, a clearly stated gap, and — importantly — the paper hands over a complete
set of validated assays and reagents to test it with (see `01_paper_summary_osterman2026.md` §4).

---

## 2. What is already established

**[FACT]** aRES binds T7 DNA polymerase directly — reciprocal His pull-downs, using
catalytically-dead aRES (R182A / E236A) to avoid killing the cells (Fig. 4H).

**[FACT]** The polymerase is *sufficient* to activate NAD⁺ cleavage: purified aRES (*K. pneumoniae*
S_34AVR) + purified T7 DNApol + NAD⁺ → ADPR-1P (Fig. 4I).

**[FACT — and this one is underappreciated]** That in vitro reaction contained **no DNA and no
dNTPs** (STAR Methods: 100 µM NAD⁺, 0.6 µM aRES, 1 µM T7 DNApol, 100 mM sodium phosphate pH 8,
1 h RT). Therefore:

> **aRES senses the polymerase *as a protein*, not its enzymatic activity, not a replication
> intermediate, and not DNA.**

That is a genuinely helpful constraint. It means a protein–protein docking approach is addressing
the real biological question, not a proxy — the system is a clean two-body recognition event.

**[FACT]** In vivo, co-expressing T7 DNApol with aRES causes growth arrest with **no phage
present** (Fig. 4G) — same conclusion, independent route.

**[FACT]** Two independent T7 escape mutants each carry a single point mutation in gp5:
**S122P** and **D106G** (Fig. 4F). Both escape *both* the *E. coli* UC5 and *K. pneumoniae*
S_34AVR aRES proteins.

**[FACT]** PF18870 is fused to HEPN, to RES, and occasionally to an ABC ATPase; some RES-fused
versions carry an extra N-terminal Zn-ribbon. It is ~100 aa and mostly α-helical.
([NCBI CDD pfam18870](https://www.ncbi.nlm.nih.gov/Structure/cdd/pfam18870);
[Anantharaman et al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3710099/) Table 2.)

---

## 3. The main structural lead: both escape mutations land in gp5 "insertion1"

This is the most useful thing to come out of the literature search, and it is **not** discussed
in the paper.

**[FACT]** T7 gp5 residues **~102–122** form a discrete structural element called **insertion1** —
a helix-loop-helix motif in the **3′→5′ proofreading exonuclease domain**. It makes stabilizing
contacts that shuttle the primer strand into the exonuclease active site.
(From the crystal structure of full-length T7 DNA polymerase with thioredoxin and dsDNA; see
`04_data_inventory.md` for the primary refs to pull — **this boundary needs confirming against the
primary structure paper, I have it second-hand**.)

**[INFERENCE] Both escape mutations fall inside this single small element:**

```
T7 gp5 exonuclease domain
├─ ... 
├─ insertion1   residues ~102 – 122   (helix–loop–helix)
│     ├─ D106G   ← escape mutation #1
│     └─ S122P   ← escape mutation #2
└─ ...
```

Two *independently isolated* escape mutations converging on one ~20-residue surface element is a
strong signal. The straightforward reading:

> **P1 — Primary prediction: insertion1 of the T7 gp5 exonuclease domain is the epitope recognized
> by aRES.**

Three things make this more than a coincidence:

1. **The escape phages are viable.** They replicate well enough to be isolated and plaque. So the
   mutations do not break polymerase function — exactly what you expect from a surface-epitope
   mutation, and *not* what you expect if the mutations acted by destabilizing the protein.
2. **Chemistry fits.** S122P inserts a proline at the end of a helix (classic helix breaker,
   distorts a defined local fold); D106G removes a charged side chain *and* adds backbone
   flexibility. Both are the sort of change that destroys a read-out surface while leaving the
   domain folded.
3. **Location fits.** These are in the exonuclease domain, far from the polymerase active site and
   far from the thioredoxin-binding domain (~258–333). An epitope, not a catalytic residue.

### Immediate practical consequence

You already have an experimentally-defined epitope on the **phage** side. That converts your blind
HADDOCK run from "hope the scoring function is right" into something you can actually
**falsify**:

> **Any docking pose that does not place aRES near gp5 residues 102–122 is probably wrong.**

Use this as a *filter on blind output*, not (initially) as a restraint — see §6.

### Honest caveats **[INFERENCE]**

- The 102–122 boundary is second-hand; confirm it from the primary structure paper before building
  restraints on it.
- Escape mutations can act indirectly: insertion1 is a proofreading element, so S122P/D106G phages
  may be mutators, and conceivably escape through some downstream effect rather than loss of
  binding. **This is directly testable** — pull-down of aRES with gp5 S122P / D106G. If the escape
  mutants still bind aRES, P1 is dead and the mutations act on activation rather than recognition.
  That experiment is cheap and should be early.
- Two mutations is a small sample. They came from a selection, so they mark *sufficient* escape
  routes, not the complete epitope.

---

## 4. ⚠️ The thioredoxin problem — check this before trusting your current docking run

**[FACT]** T7 DNA polymerase is a **1:1 heterodimer**: phage gp5 (704 aa) + **host *E. coli*
thioredoxin (TrxA, 109 aa)**. TrxA is the processivity factor; it binds a ~70-residue loop (the
thioredoxin-binding domain, TBD) inserted in the thumb subdomain, roughly residues 258–333.

**[FACT]** The enzyme used in the paper's in vitro activation experiment was **NEB M0274S**, which
is supplied as the **gp5/thioredoxin complex**, not gp5 alone.

**[INFERENCE]** So the species proven to activate aRES is the **holoenzyme**. If you are docking
gp5 alone, you may be docking a protein that:

- has a **disordered / differently-positioned TBD** (the TBD loop is poorly ordered without TrxA),
- is not the form the defense system evolved to detect.

**Recommendation.** Do not drop the gp5-only run — but add gp5+TrxA as a parallel condition and
compare. Three outcomes, all informative:

| Result | Interpretation |
|---|---|
| Same interface with and without TrxA | TrxA irrelevant to recognition; gp5-only docking is fine |
| Interface involves TrxA or the TBD | aRES may sense the *assembled holoenzyme* — a smart thing for a defense system to do, since gp5 alone is inert |
| Interface only converges with TrxA present | gp5-only docking was the wrong target all along |

**[SPECULATION]** The third possibility is biologically attractive: sensing the gp5·TrxA complex
would let the cell detect a *committed, processive, host-factor-hijacking* replication machine
rather than a harmless unassembled phage protein. Defense systems detecting "hijacking of a host
component" is a recurring theme. Worth keeping in mind, not worth assuming.

---

## 5. Why "modular sensor" (H2) is the interesting version

**[FACT]** PF18870 appears fused to at least three unrelated effectors: RES (NADase), HEPN
(RNase), ABC ATPase.

**[INFERENCE]** That distribution is the signature of a sensor module that has been **shuffled
onto different warheads** by recombination — the same architectural logic as other anti-phage
systems where a detection module and an effector module evolve semi-independently.

This yields a sharp, high-value prediction:

> **P2 — Domain-swap prediction. Chimeric aRES proteins in which NTD-1 is exchanged between
> homologs will swap phage-range specificity along with NTD-1.**

Concretely: take aRES_UC5 (*E. coli*) and aRES_GN04826 (*E. chengduensis*), which differ in phage
range (Fig. 1B heatmap). Build NTD1(UC5)–RES(GN04826) and NTD1(GN04826)–RES(UC5). Plaque-assay
both against the Fig. 1B phage panel.

- Specificity follows **NTD-1** → H1 and H2 both strongly supported. This is the clean result.
- Specificity follows **RES** → H1 is wrong; recognition lives in the C-terminal domain.
- Neither chimera works → inconclusive (folding failure); fall back on point mutants.

**This is, in my view, the single best experiment in the whole project.** It is a genetics
experiment with a binary readout, it needs only cloning + the plaque assay the lab already runs,
it does not depend on any structure prediction being correct, and it tests the central claim
directly. Everything computational below should be read as *support for, and prioritization of,
this experiment* — not as a replacement for it.

**[SPECULATION]** If H2 holds, the HEPN-NTD1 proteins are an unexplored family of phage-triggered
RNase defense systems, activated by phage proteins in the same way. That is a whole second paper,
and it is the reason to care about this domain beyond aRES.

---

## 6. On the blind HADDOCK run currently going

Straight assessment, since you asked for one:

**A truly blind, ab-initio HADDOCK run on aRES (459 aa) × gp5 (704 aa) is unlikely to give a
trustworthy answer on its own.** That is ~1,160 residues with no restraints; the sampling problem
is severe and HADDOCK's scoring function is not built for ab-initio search at that scale. Expect
poorly separated clusters.

That does **not** mean kill the run. It means:

1. **Keep it as the unbiased control.** If the blind top clusters land on insertion1 anyway, that is
   independent confirmation and far more convincing than any restrained run could ever be. Note
   this now, before you look — it is only a prediction if you commit to it first.
2. **Score it against P1 regardless of rank.** Don't only look at cluster 1. Ask: across all
   clusters, how is contact frequency with gp5 residues 102–122 distributed vs. the rest of the
   surface? And on the aRES side, NTD-1 vs. RES domain? A contact-frequency map over all models is
   far more informative than the top pose. Script stub: `analysis/scripts/` (to be written once
   structures are in).
3. **Run AlphaFold3 co-folding in parallel — this should be your primary method, not HADDOCK.**
   The paper's own successful approach for aRES–aRad1 was AF3 co-folding, reaching ipTM 0.88,
   confirmed experimentally. AF3 is free on the [AlphaFold Server](https://alphafoldserver.com/),
   handles this size, and you have a published positive control to validate against (§T0 of the
   plan). Use HADDOCK afterwards for **refinement** of an AF3 pose, or for the
   *informed* run with insertion1 restraints — which is what HADDOCK is genuinely good at.

See `03_computational_plan.md` for the ordered plan.

---

## 7. Mechanism of activation — open, and barely addressed

**[FACT]** Nothing in the paper explains *how* binding activates the NADase.

Candidate mechanisms **[SPECULATION]**, listed with what would distinguish them:

| Model | Prediction that distinguishes it |
|---|---|
| **Autoinhibition release.** NTD-1 occludes the RES active site; polymerase binding pulls it away. | Isolated RES domain (no NTD-1) should be constitutively toxic. AF3 model of full-length aRES should show NTD-1 packed against the RES pocket. |
| **Allosteric relay.** Binding at NTD-1 reorganizes the RES active site at a distance. | Full-length apo model shows NTD-1 *away* from the pocket; conserved NTD-1↔RES interface residues not at the active site. |
| **Oligomerization.** Binding drives aRES dimer/oligomer assembly that completes the active site. | SEC-MALS / crosslinking ± polymerase. AF3 of aRES homodimer. |
| **Substrate/cofactor gating.** NAD⁺ access is sterically blocked until binding. | Compare NAD⁺-bound AF3 model (ipTM 0.95, Fig. 4C) pocket accessibility in the context of full-length aRES. |

**The autoinhibition model is the cheapest to test and should go first**, because the test is a
single construct: express the RES domain alone and ask whether it is toxic without any phage.
If it is → autoinhibition, and NTD-1's job is both sensing *and* restraining. That would be a
satisfying, coherent answer and it costs one cloning step.

Note this is also a **confound for the chimera experiment (P2)**: if NTD-1 and RES co-evolved an
autoinhibitory interface, cross-species chimeras may be constitutively toxic or dead rather than
specificity-swapped. Build the chimeras, but include the RES-domain-alone construct so you can
interpret a toxic chimera correctly.

---

## 8. Alternative hypotheses that must be excluded

Listed so they are on record before the data comes in.

1. **Recognition lives in the RES domain, not NTD-1.** The aRad1 inhibitor binds the *RES* active
   site; it is not impossible that the polymerase does too, and that aRad1 works by direct
   competition with the activator. → Tested by P2 (chimeras), by the NTD-1-only vs. RES-only
   co-folding runs, and by asking whether aRad1 and T7 DNApol compete in a pull-down.
2. **Both domains contribute** — a composite interface spanning NTD-1 and RES. Quite likely, and
   would still support the spirit of H1. The isolated-domain docking runs distinguish this.
3. **The escape mutations act on activation, not binding.** → Pull-down with gp5 S122P/D106G (§3).
4. **aRES recognizes something else and gp5 is one of several triggers.** ✅ **This is now
   established, not merely possible.** λ-vir encodes **no DNA polymerase at all**, yet UC5 and
   GN04826 each give ~1,000-fold defense against it. So at least one aRES trigger is not a
   polymerase, and "aRES senses the phage DNA polymerase" holds for T7 but is not the general rule.
   → **`06_phage_range_matrix.md` §2**, which proposes that NTD-1 may instead sense phage proteins
   that *hijack host replication factors* (T7 gp5 captures thioredoxin; λ P captures DnaB), and
   gives the co-expression experiment that tests it.
5. **Indirect recognition** via a host factor. Argued against by the DNA-free, two-component in
   vitro reconstitution **[FACT]**, but TrxA is a host factor and is present in the NEB enzyme —
   see §4. Not fully excluded.

---

## 9. Falsifiable predictions, collected

| # | Prediction | Test | Cost |
|---|---|---|---|
| **P1** | aRES contacts gp5 insertion1 (~102–122) | AF3 co-folding; blind HADDOCK contact map | low |
| **P2** | NTD-1 swap transfers phage specificity | chimera + plaque assay | **medium — highest value** |
| **P3** | gp5 S122P/D106G lose aRES binding in pull-down | pull-down, as Fig. 4H | low |
| **P4** | Conserved surface patch on NTD-1, away from RES active site, = binding site | conservation mapping | low (compute only) |
| **P5** | Mutating that patch kills defense but preserves NADase | plaque assay + NAD⁺ assay, vs. R182A control | medium |
| **P6** | RES domain alone is constitutively toxic (autoinhibition) | single construct, toxicity assay | low |
| **P7** | HEPN-NTD1 proteins are also phage-triggered defense systems | bioinformatics first, then cloning | high |
| **P8** | gp5 **S122A** discriminates the serine hydroxyl from the proline's backbone rigidity as the reason S122P escapes | co-expression toxicity (Fig. 4G) + pull-down (Fig. 4H) | **low — best ratio in the project** |
| ~~P9~~ | *withdrawn.* Tested a single cofold-derived residue pair by mutating the aRES side of it. The pair was not reproduced and is dropped ([`09` §4](09_blind_zn_run_results.md)); the number is left vacant rather than reused | — | — |
| **P10** | aRES **E89A / E93A / D90A** and gp5 **R111A** lose responsiveness while retaining NADase | 4G/4H + NAD⁺ assay, vs R182A/E236A controls | low |

**P8 needs no model.** S122P is one of the paper's two escape mutations, so asking whether the
hydroxyl or the proline's backbone rigidity is responsible is a question the genetics already poses.
See `05_status_after_vault_integration.md` §4 for why it needs no mutant phages.

**P10 is the first functional test ever aimed at PF18870**, and it tests the charge-complementarity
model — acidic NTD-1 against the loop's basic patch — rather than one predicted contact. E89 is
acidic in all three orthologs, E93 loses its charge in both, and D90 is invariant, so the three ask
different questions and are not interchangeable
([`09` §5](09_blind_zn_run_results.md), [`10` §6](10_ortholog_zn_cofolds.md)).

**P5 carries a critical control design note:** an NTD-1 mutant that abolishes defense is only
interesting if the protein is still a competent NADase. You must show loss of *defense* with
retained *catalysis* — otherwise you have just made another folding mutant. Use the paper's
R182A/E236A as the catalysis-dead reference points, and the in vitro NAD⁺ assay (Fig. 4I format,
without polymerase it should be inactive; the right positive control is a lysate-activation or a
constitutively-active construct).

---

## Sources

- [Osterman et al. 2026, *Cell Host & Microbe* 34:991–999](https://doi.org/10.1016/j.chom.2026.05.004)
- [NCBI CDD pfam18870 (HEPN_RES_NTD1)](https://www.ncbi.nlm.nih.gov/Structure/cdd/pfam18870)
- [Anantharaman, Makarova, Burroughs, Koonin, Aravind 2013, *Biology Direct* 8:15](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3710099/)
- [T7 DNA polymerase — overview](https://en.wikipedia.org/wiki/T7_DNA_polymerase)
- [Richardson lab publications (T7 replication structures)](https://richardson.med.harvard.edu/publications_190-231/230.pdf)
