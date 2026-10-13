# Osterman et al. 2026 — aRES: structured summary

**Bacterial defense via RES-mediated NAD⁺ depletion is countered by phage phosphatases**
Osterman I, Hurieva B, Yirmiya E, Leavitt A, Hochhauser D, Itkin M, Malitsky S, Sorek R.
*Cell Host & Microbe* **34**, 991–999, June 10 2026. DOI: [10.1016/j.chom.2026.05.004](https://doi.org/10.1016/j.chom.2026.05.004)

Sorek lab (Weizmann). Lead contact: Rotem Sorek. First author / correspondence: Ilya Osterman.

> Scope note: this summary covers the main text + STAR Methods of the PDF. Supplementary
> figures (S1–S8) and tables (S1–S4) were **not** in the provided file and still need to be
> obtained — see `04_data_inventory.md`. Several of the most important items for this project
> (the aRES homolog list, the Fig. 1B phage-range values) live in those supplements.

---

## 1. The system in one paragraph

aRES ("anti-phage RES") is a **standalone, single-protein** bacterial defense system built from two
Pfam domains. On phage infection it is **activated by direct binding to a phage protein** and then
degrades cellular NAD⁺, causing abortive infection. Crucially, aRES cleaves NAD⁺ in a
**non-canonical** way: it produces **nicotinamide (Nam) + ADPR-1″-phosphate (ADPR-1P)** rather than
Nam + ADPR. Because ADPR-1P is not a substrate for the phage NAD⁺-reconstitution enzyme Adps,
the common phage counter-defense pathway NARP1 cannot rescue the phage. Phages fight back with
(a) an extended NARP1 operon carrying a Macro-domain phosphatase (Adpp) that converts ADPR-1P → ADPR,
and (b) aRad1, a small protein that binds and blocks the aRES active site.

## 2. Domain architecture — the centre of this project

aRES (*E. coli* UC5) is **459 aa**, with two annotated Pfam domains (Fig. 1B):

| Region | Pfam | Name | Status |
|---|---|---|---|
| N-terminal | **PF18870** | "HEPN/RES N-terminal domain 1" (NTD-1) | **uncharacterized** |
| C-terminal | **PF08808** | RES | catalytic NADase (Fig. 1B; the Methods text writes PF00808 — see `docs/08`) |

The authors' own closing statement is the opening for this project:

> "Presumably, the N-terminal domain of the aRES protein (Pfam accession PF18870), which is
> annotated as the 'HEPN/RES N-terminal domain' due to its presence in RES-containing proteins
> (Figure 1B), is responsible for recognition of the T7 DNA polymerase."

This is stated as a **presumption, not a result.** No structure, no mutation, no binding data
localizes recognition to NTD-1 anywhere in the paper. That is the gap.

They immediately add the second open question:

> "We note, however, that some aRES proteins defend against other phages including T2, T4, and T5,
> whose DNA polymerases are dissimilar to the T7 polymerase. Future studies will be necessary to
> understand the structural and molecular basis for DNA polymerase recognition by aRES and how this
> recognition activates NAD⁺ cleavage."

## 3. Key results, figure by figure

### Fig. 1 — aRES protects via NAD⁺ depletion
- Found in **defense islands**, standalone (not operonic with an antitoxin) — argues it is *not*
  a toxin-antitoxin module (1A). Representative loci from *E. coli* UC5, *K. pneumoniae* S_34AVR,
  *E. chengduensis* GN04826, *Pseudomonas* sp. LFM046, *A. oleivorans* DQS4.
- Three homologs cloned and expressed in *E. coli* MG1655: **UC5**, **S_34AVR**, **GN04826**.
  All three protect against diverse phages (1B heatmap: T2, T4, T5, T7, SECphi4, SECphi18,
  SECphi17, SECphi27, λ-vir; fold defense up to >100,000).
  → **The heatmap is the key evidence that different aRES have different phage ranges.**
  Exact per-cell values must be read from the figure/supplement, not from this summary.
- Defense is associated with **cell death** (abortive infection), not phage-genome clearance (1C, S2A).
- NAD⁺ is depleted in aRES cells during T7 infection (1D).
- Active-site point mutants abolish defense: **R182A, R184A, R215A, V293A, E236A, S234A** (1E).
  These were chosen by structural alignment to RES toxins.
  → **R182A and E236A are the catalytically-dead-but-folded variants used for pull-downs.**

### Fig. 2 — The chemistry, and why NARP1 fails
- aRES-expressing infected cells accumulate **ADPR-1P**, not ADPR (2A).
- LC-MS of lysate reaction: NAD⁺ → **Nam + ADPR-1P** (2B, 2C). Contrast with SIR2/TIR/SEFIR, which
  give Nam + ADPR.
- NARP1 cannot cancel aRES defense; **NARP2 can** (2D). Mechanistic reason: Adps (the first NARP1
  enzyme) cannot use ADPR-1P as substrate.

### Fig. 3 — The phage counter-counter-defense (NARP1-AP)
- ~15% of NARP1-carrying phages also encode a **Macro-domain** protein, named **Adpp**.
  HHpred's best hit with known function: yeast **Poa1p**, a known ADPR-1P phosphatase.
- The three-gene operon Adps–Adpp–Namat from *Vibrio* phage Gary overrides aRES defense (3B, 3C).
  Inactivating **any** of the three abolishes anti-defense → all enzymes required.
- Purified Adpp converts ADPR-1P → ADPR in vitro (3E).
- NARP1-AP also cancels Nezha defense (Fig. S6).

### Fig. 4 — **The figure that matters most here**

**4A–4E: aRad1, a phage inhibitor of aRES**
- Discovered by large-scale **AlphaFold3-Multimer co-folding screen**: ~133,000 short phage proteins
  of unknown function from IMG/VR v4 (~112 M proteins) co-folded against the RES domain.
- Hit: a family whose members bind aRES with high confidence via a **flexible loop entering the
  aRES active-site pocket**. Seven tested, **six abolished aRES defense** → named **aRad1**.
- AF3 model aRES–aRad1b: **ipTM 0.88** (4A). The loop penetrates the RES active site; a conserved
  **W58** from another part of aRad1 overlaps the RES binding pocket; loop residues **G15, Q16**
  contact aRES **R182, R184** (4B).
- AF3 model aRES + NAD⁺: **ipTM 0.95** (4C). Same pocket, same arginines → aRad1 is an active-site
  competitor.
- Mutations abolishing aRad1 anti-defense: **W58A** and **ΔG15** (4D).
- **Pull-down:** His-aRad1 co-purifies aRES → direct physical complex. The flexible-loop mutant
  (ΔG15) pulls down aRES only weakly → the loop drives the interaction (4E).

**4F–4I: T7 DNA polymerase is the activator — the trigger for this project**
- Two independently isolated **T7 escape mutants** were whole-genome sequenced. Each carried a
  **single point mutation in the phage DNA polymerase (gp5)**:
  - **S122P**
  - **D106G**
- Both escape **both** the *E. coli* UC5 and the *K. pneumoniae* S_34AVR aRES (4F)
  → **the recognition determinant is shared between these two divergent aRES proteins.**
- **Co-expression of T7 DNApol with aRES causes growth arrest with no phage present** (4G) — the
  polymerase alone is sufficient to trigger toxicity. The S122P polymerase does not (or much less).
- **Pull-down:** His-tagged aRES (catalytically dead R182A or E236A, to avoid killing the cells)
  pulls down T7 DNA polymerase, and reciprocally (4H) → **direct physical interaction.**
- **In vitro reconstitution (4I):** purified *K. pneumoniae* S_34AVR aRES + purified T7 DNA
  polymerase (NEB M0274S) + NAD⁺, 1 h, RT → **ADPR-1P produced**. Without the polymerase, no
  turnover. → The polymerase alone is **sufficient** for activation; no other phage or host factor
  needed.

## 4. Methods details that constrain the experimental plan

Taken from STAR Methods — these are practical facts worth knowing *before* the rotation starts.

- **aRES from *E. coli* UC5 could not be purified as a soluble protein.** The authors switched to
  the *K. pneumoniae* S_34AVR homolog (C-terminal 6×His, pBAD, LOBSTR-BL21(DE3)-RIL, induction
  0.5 mM IPTG overnight at 16 °C). **→ For any in vitro work, start with S_34AVR, not UC5.**
- Adpp was expressed as an N-terminal His14-bdSUMO fusion (pET28-His-bdSUMO), cleaved with bdSENP1.
- **T7 DNA polymerase used in vitro was the commercial NEB M0274S.** That product is the
  **gp5 + *E. coli* thioredoxin complex**, not gp5 alone — see `02_hypothesis_NTD1_sensor.md` §4.
- Pull-downs: Ni-NTA magnetic beads, elution in 500 mM imidazole, SDS-PAGE readout. aRES/T7 DNApol
  co-expression used 1 mM IPTG + 0.2% arabinose, 25 °C, mid-log.
- Growth/infection assays at **25 °C** in MMB (LB + 0.1 mM MnCl₂ + 5 mM MgCl₂), OD₆₀₀ every 10 min,
  Tecan Infinite 200.
- NAD⁺ quantification: Promega NAD/NADH-Glo.
- Structure prediction: **AlphaFold3 with default parameters**, via the public
  [AlphaFold Server](https://alphafoldserver.com/). NAD⁺ was supplied as SMILES (string given in
  Methods, reproduced in `data/external/README.md` once available).

### The authors' own AF3-Multimer screening protocol (reusable here)

This is the protocol that found aRad1, and the template for any screen in this project:

1. Phage protein structures predicted with **AlphaFold2 v2.1.1**; clusters from **MMseqs2** with
   ≥30 non-identical members; representative per cluster; length filter **300–600 aa**; proteins
   with HHpred hits to known function discarded. Yield: 182,179 structures.
2. **AF3** co-folding of each phage sequence with the RES domain.
3. Top 10 by ipTM inspected manually; the one with a flexible loop in the active-site pocket taken
   forward.
4. **Foldseek** (release 5.53465f065) against the 182,179 set, probability > 90%, expanded to full
   MMseqs2 clusters → 156 homologs.
5. Each of the 156 co-folded with **five different RES domains** (IMG accessions
   **2508547945, 2574243703, 2703310322, 2729545185, 2753914621**), **5 AF3 models each**,
   interaction score = **mean ipTM across the 5 models**.
6. Seven top homologs taken to experiment.

**Note the asymmetry:** the screen co-folded phage proteins against the **RES domain**, i.e. it was
designed to find active-site *inhibitors*. A screen against **NTD-1** has, as far as this paper
shows, **never been run.** That is a concrete, cheap, novel experiment — see `03_computational_plan.md` §T5.

## 5. Detection of aRES systems (for homolog collection)

- Source: ~135 M genes from ~38,000 bacterial/archaeal genomes (IMG, Oct 2017).
- Candidates: annotated **PF00808** *(sic — the Methods text says PF00808 throughout, while Fig. 1B
  labels the domain PF08808. PF08808 is RES; PF00808 is a histone-like family. Quoted here as the
  paper writes it; see `docs/08` for which to actually use.)*, or in protein clusters where that
  annotation was the most abundant; length filter 300–600 aa.
- MSAs via hhblits/UniRef30 (`-oa3m -n 3 -M first`), searched with hhsearch against Pfam v35
  (`-p 80 -loc -z 1 -b 1 -ssm 2 -sc 1 -seq 1`).
- **Final criterion: hits to BOTH PF00808 *(as written; read PF08808)* and PF18870 at p > 80%.**
- Result: aRES homologs in **~4% of genomes** — 1,526 of 38,149 (Fig. S7A, Tables S2 and S3).

→ **Tables S2/S3 are the homolog list this project needs for conservation analysis.**

## 6. Numbers worth remembering

| Quantity | Value |
|---|---|
| aRES (UC5) length | 459 aa |
| aRES prevalence | ~4% of 38,149 genomes (1,526) |
| aRad1 ipTM (AF3) | 0.88 |
| aRES + NAD⁺ ipTM (AF3) | 0.95 |
| Fold defense (plaque assay) | up to >100,000 |
| T7 escape mutations in gp5 | S122P, D106G |
| aRES catalytic/defense-null mutants | R182A, R184A, R215A, V293A, E236A, S234A |
| aRad1 anti-defense-null mutants | W58A, ΔG15 |
| NARP1-carrying phages with Macro domain | ~15% |
| Bacteria defending via NAD⁺ depletion | now estimated >20% |

## 7. What the paper does *not* establish (the project's room to work)

1. **No evidence that NTD-1 does the recognizing.** It is a presumption from domain annotation.
2. **No structure of the aRES–T7 DNApol complex** — predicted or experimental. Only pull-downs.
3. **No mapping of the aRES surface** that contacts the polymerase. No aRES interface mutants.
4. **No mechanism for activation.** Binding → NADase activation is a black box. Is it allosteric
   release of an autoinhibited state? Does NTD-1 occlude the RES active site until displaced?
   Does binding drive oligomerization? None of this is addressed.
5. **No explanation for the T2/T4/T5 phenotype** — different aRES defend against phages with
   dissimilar polymerases. Either NTD-1 recognizes a shared structural feature, or there are
   *other triggers* entirely.
6. **No explanation of the specificity differences** between the three tested aRES homologs.
7. **gp5 escape mutations are not structurally interpreted** anywhere in the paper.

Point 7 is where the first real lead came from — see `02_hypothesis_NTD1_sensor.md`.

---

## Sources

- [Osterman et al. 2026, Cell Host & Microbe 34:991–999](https://doi.org/10.1016/j.chom.2026.05.004) (the PDF)
- [Pfam PF18870 / NCBI CDD entry](https://www.ncbi.nlm.nih.gov/Structure/cdd/pfam18870)
- [Anantharaman et al. 2013, Biology Direct 8:15 — HEPN superfamily](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3710099/)
