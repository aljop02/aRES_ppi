# Computational plan

Ordered by **information gained per unit of effort**. Tiers T0–T2 are the core and could realistically
be done in the first week or two; T3–T6 depend on what they show.

Principle throughout: **every prediction method must first reproduce a known answer.** The paper
handed us two published positive controls with numbers attached (aRES–aRad1, ipTM 0.88; aRES–NAD⁺,
ipTM 0.95). Use them. A pipeline that cannot reproduce those cannot be trusted on the unknown case.

---

## T0 — Validate the pipeline on the published positive controls

**Do this before anything else. It is the difference between a result and a guess.**

| Run | Components | Expected (from paper) |
|---|---|---|
| T0.1 | aRES (*E. coli* UC5, full length) + aRad1b | **ipTM ≈ 0.88**; aRad1 flexible loop in the RES active-site pocket; aRad1 **G15/Q16** contacting aRES **R182/R184**; conserved **W58** overlapping the pocket |
| T0.2 | aRES (UC5) + NAD⁺ (SMILES from STAR Methods) | **ipTM ≈ 0.95**; NAD⁺ in the RES pocket, H-bonds to R182/R184 |
| T0.3 | aRES (UC5) + aRad1b **ΔG15** | weaker/altered interface (Fig. 4E: pulls down aRES only weakly) |

Method: AlphaFold3, default parameters, [AlphaFold Server](https://alphafoldserver.com/) — exactly
what the authors used. 5 models per run; report mean and spread of ipTM, not a single number.

**Acceptance criterion:** T0.1 within ~0.1 ipTM of 0.88 **and** the contact residues above
recovered. If you cannot reproduce this, stop and debug (wrong sequence? wrong aRad1 homolog?
chain stoichiometry?) before proceeding. Record everything in `analysis/afold/T0_validation/`.

> Why this matters more than it sounds: it also calibrates what ipTM value means *for this
> particular complex type*. If your T1 run gives 0.6, you need to know whether 0.6 is "weak" in a
> system where the known true positive scores 0.88.

---

## T1 — Predict the aRES–T7 DNA polymerase complex

The central run. Use **S_34AVR (*K. pneumoniae*) as the primary aRES** — that is the homolog with
the direct in vitro evidence of activation by T7 DNApol (Fig. 4I) *and* the one that is actually
purifiable (UC5 is not — see `01_paper_summary_osterman2026.md` §4).

| Run | Chain A | Chain B | Purpose |
|---|---|---|---|
| T1.1 | aRES **S_34AVR** full | T7 gp5 full | primary prediction, best ground truth |
| T1.2 | aRES **UC5** full | T7 gp5 full | the paper's main strain; cross-check |
| T1.3 | aRES S_34AVR full | T7 gp5 **+ TrxA** | holoenzyme — the species actually used in vitro (§4 of hypothesis doc) |
| T1.4 | aRES **NTD-1 only** | T7 gp5 full | **direct test of H1** |
| T1.5 | aRES **RES domain only** | T7 gp5 full | control; should be *worse* if H1 holds |
| T1.6 | aRES S_34AVR full | T7 gp5 **exonuclease domain only** (~1–180) | tests P1 with a reduced search space |

5 AF3 models each.

**Domain boundaries for T1.4/T1.5 — use the PAE-derived ones, not the Pfam ones.** From the AF3 PAE
analysis in `notes/vault/` (see `05_status_after_vault_integration.md` §1.1):

| Construct | Residues |
|---|---|
| **NTD-1 / N-block** | **1–120** (or 1–124) |
| **RES side** | **126–447** (retains R182, R184) |
| Full length | 1–459 |
| Tail-trimmed | 1–447 |

InterPro calls PF18870 at **49–136** and PF00808 at **199–363**, but the structural module runs
1–120 with a hinge at 120–125. Cutting at the Pfam coordinates would slice through a folded block
at one end and through the hinge at the other. Numbering is for ***E. coli* UC5**.

### What to record for each run

Not just ipTM. Build a per-run row with:

- mean + SD ipTM across 5 models, pTM
- **interface PAE** (the mean PAE of inter-chain residue pairs within 8 Å) — more diagnostic than
  ipTM for localizing *which part* of the interface is confident
- contact list (heavy-atom < 5 Å), split by aRES domain: NTD-1 vs. RES vs. linker
- **fraction of gp5 interface residues falling in 102–122** ← the P1 readout
- whether the RES active site (R182/R184/E236 numbering is UC5 — map to S_34AVR!) is occluded

> ⚠️ **Residue numbering trap.** All the published mutation numbers (R182, R184, E236, S234, R215,
> V293) are for ***E. coli* UC5 aRES**. S_34AVR is a different protein with different numbering.
> Build a UC5↔S_34AVR alignment **first** and keep a mapping table in
> `analysis/conservation/numbering_map.tsv`. Getting this wrong will silently corrupt every
> downstream interpretation.

### Decision point

- Interface is predominantly **NTD-1**, and lands on gp5 insertion1 → **H1 + P1 supported.** Go to T2.
- Interface is predominantly **RES domain** → H1 challenged; the alternative in hypothesis doc §8.1
  is live. Re-plan; the chimera experiment (P2) becomes even more important as the arbiter.
- **No convergent interface** across models/seeds → AF3 cannot do this one. Do not force it. Go to
  T4 (conservation) — which is independent of structure prediction — and lean on the wet-lab
  chimera experiment.

A null result here is informative and should be written up as such, not buried.

---

## T2 — The escape-mutant discrimination test

Re-run T1.1 / T1.2 / T1.3 with **gp5 S122P** and **gp5 D106G**.

**Prediction if P1 holds:** interface degrades — lower ipTM, worse interface PAE, or the pose moves
off insertion1.

### ⚠️ Read this before interpreting the result

**AlphaFold is notoriously insensitive to single point mutations.** It frequently returns a nearly
identical structure and score for a mutant that is biologically dead. Therefore:

- A **positive** result (clear degradation) is meaningful and encouraging.
- A **null** result (no change) is **not** evidence against P1. Do not report it as such.

To get more signal than ipTM alone:

1. Compare **interface PAE maps**, not just scalar scores — local confidence around 102–122 is more
   sensitive than a global score.
2. Compare **contact-frequency maps across all 5 models** rather than the top model.
3. Consider an orthogonal estimate: MD-based or FoldX/Rosetta ΔΔG of binding on the T1 pose. These
   have their own large error bars, but they fail differently from AF3, and agreement between two
   methods that fail differently is worth something.
4. The real test of this is the wet-lab pull-down (P3), which is cheap. Treat T2 as a prioritization
   exercise, not a verdict.

---

## T3 — HADDOCK, done properly

Three runs, in increasing order of information injected. Keep all three — the comparison *between*
them is itself the result.

**T3.1 — the blind run (already going).** Keep it as the unbiased control. Analyze it as described
in hypothesis doc §6: contact-frequency map over *all* clusters, not just the top pose.

**T3.2 — epitope-restrained.** Ambiguous Interaction Restraints (AIRs):
- **gp5 active residues:** insertion1 (~102–122), filtered to those with relative SASA > 40% in the
  monomer structure.
- **gp5 passive:** surface neighbours within ~6.5 Å of the active set (HADDOCK can auto-define).
- **aRES active:** from T4 conservation analysis (conserved, surface-exposed NTD-1 patch). If T4 is
  not done yet, run with gp5-side restraints only — a one-sided AIR setup is legitimate and is the
  standard approach when you know one epitope.

**T3.3 — refinement of the AF3 pose.** Take the best T1 model, run HADDOCK in refinement mode
(`refine` / water refinement). This is where HADDOCK genuinely adds value over AF3: side-chain
repacking and a physics-based energy, which AF3 does not give you.

**Input structure note.** Use the real experimental structures where they exist rather than
predictions — for T7 DNA polymerase there are crystal structures of the gp5·TrxA·DNA complex
(see `04_data_inventory.md`). Strip the DNA. For aRES there is no structure; use the AF model, and
**check the pLDDT** — low-confidence loops dock badly and should be flagged, possibly treated as
fully flexible in HADDOCK.

---

## T4 — Conservation and specificity-determining positions

**Highest value per CPU-hour in the whole plan, and it does not depend on any docking being right.**
If one thing gets done, make it this.

### T4.1 Build the homolog set
Source: paper Tables **S2/S3** (1,526 genomes with aRES). If those tables are not obtainable,
reproduce the authors' detection (Methods §"Detection of putative aRES systems"): PF00808 candidates,
300–600 aa, hhblits/UniRef30 MSA, hhsearch vs Pfam v35, require **both** PF00808 and PF18870 at
p > 80%.

### T4.2 Two alignments, deliberately
1. **NTD-1 from aRES proteins** (PF18870 + PF00808 architecture)
2. **NTD-1 from HEPN proteins** (PF18870 + HEPN architecture)

Aligning both and comparing is what tests **H2**. If NTD-1 is a modular sensor, the two sets should
share a conserved *surface* (the sensing face) while differing where they pack against their
respective effectors. **That comparison is the single cleanest bioinformatic test of H2** and I have
not seen it done anywhere.

MAFFT L-INS-i for ≤1,000 seqs; HMMER for profile building; trim with trimAl, but keep an untrimmed
copy — trimming can delete exactly the variable specificity positions you are looking for.

### T4.3 Score and map
- Per-position conservation: Rate4Site / ConSurf, or Jensen–Shannon divergence.
- Map onto the AF model of aRES NTD-1; compute per-residue relative SASA.
- **Target set = conserved ∧ surface-exposed ∧ not in the NTD-1↔RES interface.** That is the
  candidate recognition surface. This is P4.

### T4.4 Specificity-determining positions (SDPs)
The point of the project is that *different aRES recognize different phage proteins*. So look for
positions that are **conserved within subfamilies but differ between them** — classic SDP analysis
(e.g. SDPfox, GroupSim, or mutual information with the subfamily partition).

Anchor the partition with the three experimentally characterized homologs (UC5, S_34AVR, GN04826)
and their phage ranges from Fig. 1B. Three anchors is thin, so treat subfamily assignment as a
clustering result to be sanity-checked, not ground truth.

**Convergence test:** SDPs ∧ the T1 predicted interface. If independent evolutionary signal and
independent structural prediction point at the same residues, that is a genuinely strong result —
much stronger than either alone, because the methods share no assumptions.

---

## T5 — The T2/T4/T5 puzzle, and a screen nobody has run

**The puzzle [FACT]:** some aRES defend against T2, T4 and T5, "whose DNA polymerases are dissimilar
to the T7 polymerase."

Two explanations, both testable:

**T5.1 — Same trigger class, divergent sequence.** Do T2/T4/T5 polymerases carry an insertion1-like
element despite low sequence identity? Structure-based search: extract the gp5 insertion1 motif,
**Foldseek** it against predicted structures of T2/T4/T5 polymerases. Sequence identity would miss
this; structural similarity would not.

**T5.2 — Different trigger entirely.** Then comparative genomics should find it. You have a
defense-range matrix in Fig. 1B (3 aRES × ~9 phages). Ask: which phage proteins are present in the
defended phages and absent or divergent in the non-defended ones? With 9 phages the statistics are
weak, but it nominates candidates for co-folding, which is the real filter.

**T5.3 — ⭐ Run the AF3 screen the authors never ran.**

This is the clearest novelty opportunity in the project. Re-read the Methods: the authors' 182,179-structure
screen co-folded phage proteins **against the RES domain** — it was designed to find active-site
*inhibitors*, and it found aRad1. **Nobody has screened phage proteomes against NTD-1.**

A screen against NTD-1 asks the opposite question: *what phage proteins does this sensor bind?*
Start small and tractable rather than attempting 182,179 structures:

- Phage proteomes of the Fig. 1B panel (T2, T4, T5, T7, SECphi4/17/18/27, λ) — a few hundred to a
  few thousand proteins. Entirely feasible.
- Co-fold each against NTD-1 from each of the three characterized aRES.
- Rank by mean ipTM over 5 models (the authors' own criterion).
- **Built-in positive control: T7 gp5 must come out near the top for an aRES that defends against
  T7.** If it does not, the screen is not working — and you know that before interpreting anything
  else. That is a well-posed screen.

Compute cost is the limiting factor; the AlphaFold Server is rate-limited, so scope this to what is
actually runnable (consider a local AF3/Boltz-style setup, or triage by size/annotation first as the
authors did with their 300–600 aa filter).

---

## T6 — Deliverable for the bench

The plan's output should be a **ranked, annotated construct list**, not a pile of PDB files. Target
format, one row per construct, in `results/`:

| Construct | Rationale | Predicted effect | Readout | Priority |
|---|---|---|---|---|

Construct classes to populate:

1. **NTD-1 surface mutants** (from T4.3 ∩ T1) — predicted to abolish defense, retain NADase.
   Must be paired with the catalysis controls R182A/E236A. See hypothesis doc §9 on why the
   "retains catalysis" half is non-negotiable.
2. **Domain-swap chimeras** (P2) — NTD-1 exchanged between UC5 / S_34AVR / GN04826. **Highest
   priority of anything in this document.** Include the RES-domain-alone construct alongside, so a
   constitutively toxic chimera can be interpreted (hypothesis doc §7).
3. **gp5 insertion1 mutants** beyond S122P/D106G — predicted additional escape routes. A clean
   prospective prediction: nominate them computationally *before* testing.
4. **gp5 S122P / D106G pull-downs** (P3) — binding vs. activation. Cheap, early, high information.
5. **Truncations** — NTD-1 alone, RES alone — for the autoinhibition question (P6).

Assays all exist in the paper and need no development: plaque assay fold-defense, liquid infection
OD₆₀₀, NAD/NADH-Glo, LC-MS for ADPR-1P, Ni-NTA pull-down, and the two-component in vitro
reconstitution.

---

## Suggested order

```
T0  validate pipeline           ──► must pass before anything else
 │
 ├─ T1  aRES × gp5 co-folding   ──┐
 │   └─ T2  escape mutants        │
 │                                ├──► T6 construct list ──► bench
 ├─ T4  conservation / SDP      ──┤     (chimeras first)
 │                                │
 └─ T3  HADDOCK (blind + AIR)   ──┘
 
 T5  broader trigger question / NTD-1 screen  ──► second phase, higher risk & reward
```

**In parallel from day one:** the chimera cloning (P2). It does not depend on any computational
result, it is the most direct test of the central hypothesis, and cloning + plaque assays have a
long lead time. Do not let it wait on the docking.
