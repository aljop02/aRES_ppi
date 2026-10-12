# Data inventory — what we need, what we have, how to get it

**Current state: we have the paper's main text and nothing else.** No sequences, no structures, no
supplements are in this repo yet. This document is the acquisition checklist.

---

## ⚠️ Network restriction in the cloud session

The Claude cloud environment's network policy **denies outbound HTTPS** to every sequence and
structure database tested:

```
rest.uniprot.org        403  connect_rejected
files.rcsb.org          403  connect_rejected
www.ebi.ac.uk           403  connect_rejected
eutils.ncbi.nlm.nih.gov 403  connect_rejected
```

Web *search* works (it runs server-side), so background literature is reachable; direct database
downloads are not.

**Two ways to fix this, pick either:**

1. **Allow the domains.** Cloud environment menu in the session title bar → *Edit* → *Network
   access*: either raise the access level, or add these under **Allowed domains** (leave "Allow
   package managers" ticked):
   ```
   rest.uniprot.org
   www.ebi.ac.uk
   files.rcsb.org
   eutils.ncbi.nlm.nih.gov
   alphafold.ebi.ac.uk
   www.ncbi.nlm.nih.gov
   ```
   Docs: https://code.claude.com/docs/en/cloud-environments#network-access
2. **Upload the files into the chat**, or commit them to this repo from your own machine. Put
   anything you upload in `data/external/`.

**Nothing in this repo is a stand-in for real data.** No sequence, structure or alignment here has
been invented; where a value is uncertain it is marked "VERIFY".

---

## 1. Protein sequences needed

### Phage side

| Protein | Accession | Status | Notes |
|---|---|---|---|
| T7 DNA polymerase, gp5 | **UniProt P00581** (DPOL_BPT7) | ✅ accession confirmed | *Enterobacteria phage T7*, taxon 10760. The docking target. |
| *E. coli* thioredoxin 1 (TrxA) | UniProt **P0AA25** — VERIFY | ⚠️ accession unverified | Processivity factor; part of the NEB enzyme. See hypothesis doc §4. |
| aRad1b | from paper / IMG | ❌ not located | Needed for the T0 pipeline validation — **blocking for T0.1**. Likely in Table S1 or S4, or derivable from the IMG phage protein ID `Ga0374944_603278_30886_31098` (read from Methods). |
| T2 / T4 / T5 DNA polymerases | from the phage genomes below | ❌ | For T5.1. |

Phage genome accessions (from the paper's Key Resources Table, clean text — reliable):

| Phage | Accession |
|---|---|
| λ-vir | NC_001416.1 |
| T2 | LC348380.1 |
| T4 | AF158101.6 |
| T5 | AY543070.1 |
| T7 | NC_001604.1 |
| SECphi4 | MT331608.1 |
| SECphi17 | LT960607.1 |
| SECphi18 | LT960609.1 |
| SECphi27 | LT961732.1 |

These nine are the Fig. 1B defense panel — the proteomes for the T5.3 screen.

### Bacterial side

| Protein | Source | Status |
|---|---|---|
| aRES *E. coli* UC5 (459 aa) | IMG scaffold `Ga0127341_109` — VERIFY | ❌ |
| aRES *K. pneumoniae* S_34AVR | IMG scaffold `Ga0618786_134` — VERIFY | ❌ **primary target** (purifiable; in vitro validated) |
| aRES *E. chengduensis* GN04826 | IMG scaffold `Ga0377618_1185` — VERIFY | ❌ |
| aRES *Pseudomonas* sp. LFM046 | IMG scaffold `Ga0077281_120` — VERIFY | ❌ |
| aRES *A. oleivorans* DQS4 | IMG scaffold `Ga0375493_11` — VERIFY | ❌ |

> **All five IMG scaffold IDs were read off the small type in Fig. 1A and may contain transcription
> errors.** Check them against the paper at full resolution before using. The five RES-domain IMG
> accessions in the Methods are clean text and more reliable:
> `2508547945`, `2574243703`, `2703310322`, `2729545185`, `2753914621`.

Easiest route to the three characterized aRES sequences is probably **Table S4** (plasmids/strains),
which should contain the cloned ORFs directly.

---

## 2. Structures needed

| Structure | ID | Status | Use |
|---|---|---|---|
| T7 DNApol·TrxA·DNA·ddNTP | **PDB 1T7P** | ✅ confirmed | The docking input. 2.2 Å, Doublié et al. 1998 *Nature* 391:251–258, PMID 9440688, DOI [10.1038/34593](https://doi.org/10.1038/34593). Contains a 698-residue T7 pol chain + 108-residue *E. coli* TrxA + 2 DNA strands. **Strip the DNA before docking.** |
| Other T7 pol structures | 1T8E, 1X9S, 1X9M, 1TK5, 1SL2 — VERIFY | ⚠️ | Alternative conformations; useful for ensemble docking. IDs seen in CATH listings, not individually checked. |
| aRES | none exists | ❌ | Must be predicted (AF3/AF2). |
| aRES–aRad1 | none deposited | ❌ | The paper's Fig. 4A model was AF3; not deposited. Must re-predict → that *is* T0.1. |

Note on 1T7P: only 698 of ~704 residues are modelled, so **expect missing termini/loops**. Check
whether insertion1 (102–122) is fully resolved — if it is not, that itself is informative (a
flexible epitope), and it changes how you should treat it in HADDOCK.

---

## 3. Supplementary material — **the highest-priority acquisition**

Not in the provided PDF. Get these from
[the paper's supplementary material](https://doi.org/10.1016/j.chom.2026.05.004):

| Item | Contents | Why it matters |
|---|---|---|
| **Table S2, S3** | the ~1,526 genomes / aRES homolog list | **Blocking for all of T4.** Without it, the conservation analysis has to be rebuilt from scratch. |
| **Table S4** | plasmids, strains, phages, primers | The actual cloned aRES and aRad1 sequences. Probably the fastest route to every sequence above. |
| **Table S1** | NARP1 operons with Macro domain | Lower priority for this project. |
| **Fig. S7A** | aRES distribution | Context. |
| **Fig. S3** | RES active-site structural alignment | How R182/R184/E236 etc. were chosen — needed for the UC5↔S_34AVR numbering map. |
| **Fig. S7D, S7E** | aRES–aRad1 interface detail | Extra validation targets for T0. |
| **Fig. 1B at full resolution** | per-phage fold-defense values | **The specificity matrix the whole project rests on.** Needed to anchor the T4.4 subfamily partition. |

If the supplements are paywalled: the Sorek lab deposits data generously and the lead contact
(Rotem Sorek) is explicitly listed for reagent/data requests. Your PI may also simply have access.

---

## 4. Reagents available commercially (for the bench half)

| Reagent | Source | Cat# |
|---|---|---|
| T7 DNA Polymerase (**gp5·TrxA complex**) | NEB | M0274S |
| NAD | Sigma | NAD100-RO 10127965001 |
| ADPR | Sigma | A0752-25MG |
| NAD/NADH-Glo | Promega | G9071 |
| pBad/His | Thermo | 43001 |

---

## 5. Software

| Tool | Use | Note |
|---|---|---|
| **AlphaFold3** | primary structure prediction | [AlphaFold Server](https://alphafoldserver.com/), free, default params — exactly what the authors used. Rate-limited. |
| HADDOCK 2.4 / 3 | docking, refinement | Web server or local |
| MAFFT, HMMER, trimAl | alignments | pip/conda installable in this container |
| Foldseek | structure search | For T5.1 |
| MMseqs2 | clustering | As in the paper's screen |
| ConSurf / Rate4Site | conservation | Server or local |
| PyMOL / ChimeraX | visualization | Local |
| Biopython | parsing, analysis | ⚠️ **not installed in this container** — `pip3 install biopython` (pypi is reachable) |

---

## 6. Acquisition checklist

Ordered by what unblocks the most downstream work.

- [ ] **Supplementary Tables S2/S3/S4** ← unblocks T4 and most sequence needs
- [ ] **Fig. 1B at full resolution** ← the specificity matrix
- [ ] aRES S_34AVR sequence ← primary docking target
- [ ] aRad1b sequence ← blocking for T0 pipeline validation
- [ ] T7 gp5 (P00581) + TrxA sequences
- [ ] PDB 1T7P
- [ ] aRES UC5 + GN04826 sequences ← for the chimera design (P2)
- [ ] Confirm insertion1 boundaries from the primary T7 structure literature
- [ ] Decide: open the network policy, or upload files

Drop uploads in `data/external/`. Keep anything downloaded from a database in `data/sequences/`
or `data/structures/` with a `SOURCE.md` recording accession, date and URL — future-you will want
to know which version of a database an alignment came from.

---

## Sources

- [Osterman et al. 2026, *Cell Host & Microbe*](https://doi.org/10.1016/j.chom.2026.05.004)
- [RCSB PDB 1T7P](https://www.rcsb.org/structure/1t7p) · [PDBj 1T7P](https://pdbj.org/mine/summary/1T7P)
- [Doublié et al. 1998, *Nature* 391:251–258](https://doi.org/10.1038/34593)
- [UniProt P00581 (GONUTS mirror)](https://gowiki.tamu.edu/wiki/index.php/BPT7:DPOL)
- [T7 DNA polymerase overview](https://en.wikipedia.org/wiki/T7_DNA_polymerase)
- [Cloud environment network access docs](https://code.claude.com/docs/en/cloud-environments#network-access)
