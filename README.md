# aRES_ppi — protein–protein interactions of the aRES anti-phage defense system

Rotation project: **how does aRES recognize phage proteins, and what activates its NADase?**

Central hypothesis: the uncharacterized N-terminal domain of aRES (**Pfam PF18870**, "HEPN/RES
N-terminal domain 1", NTD-1) is a **modular phage-protein sensor** — it binds the trigger protein
(T7 DNA polymerase gp5 in the characterized case) and thereby activates the C-terminal RES NADase.

Starting point: **Osterman et al. 2026, *Cell Host & Microbe* 34:991–999**
([10.1016/j.chom.2026.05.004](https://doi.org/10.1016/j.chom.2026.05.004)), which establishes aRES,
shows it is activated by direct binding to T7 DNA polymerase, and explicitly leaves the recognition
mechanism open.

---

## Read in this order

| Doc | What it is |
|---|---|
| [`docs/01_paper_summary_osterman2026.md`](docs/01_paper_summary_osterman2026.md) | Structured summary of the paper: results, methods constraints, and a list of what it does *not* establish |
| [`docs/02_hypothesis_NTD1_sensor.md`](docs/02_hypothesis_NTD1_sensor.md) | The hypothesis, the evidence for it, the structural lead, alternatives to exclude, numbered predictions P1–P7 |
| [`docs/03_computational_plan.md`](docs/03_computational_plan.md) | Tiered plan T0–T6, ordered by information per unit effort |
| [`docs/04_data_inventory.md`](docs/04_data_inventory.md) | Every sequence/structure/supplement needed, with accessions and acquisition status |

## The three things worth knowing up front

**1. Both T7 escape mutations land in one small structural element.**
The escape mutants are gp5 **S122P** and **D106G**. Residues ~102–122 of gp5 form **insertion1**, a
helix-loop-helix in the proofreading exonuclease domain. Two independently isolated escape mutations
converging on one ~20-residue surface element is a strong lead, and it is not discussed in the paper.
→ **P1: insertion1 is the epitope aRES reads.** Use it to filter docking output.
Details and caveats in [hypothesis §3](docs/02_hypothesis_NTD1_sensor.md).

**2. T7 DNA polymerase is a heterodimer, and the activating species was the holoenzyme.**
gp5 (704 aa) binds host *E. coli* **thioredoxin** as its processivity factor. The enzyme used in the
paper's in vitro activation experiment (NEB M0274S) is the **gp5·TrxA complex**, not gp5 alone.
Docking gp5 by itself may be docking the wrong species.
→ [hypothesis §4](docs/02_hypothesis_NTD1_sensor.md).

**3. The best experiment is genetics, not docking.**
Swap NTD-1 between aRES homologs with different phage ranges and ask whether specificity follows
NTD-1. Binary readout, uses only assays the lab already runs, and does not depend on any prediction
being correct. **Start the cloning in parallel with the computational work, not after it.**
→ P2 in [hypothesis §5](docs/02_hypothesis_NTD1_sensor.md).

## Status

- [x] Paper read and summarized
- [x] Hypothesis formalized with falsifiable predictions
- [x] Computational plan drafted
- [x] Interface-analysis tooling written and tested
- [ ] **Sequences and structures acquired** ← blocking everything below
- [ ] T0: pipeline validated against the published aRES–aRad1 control
- [ ] T1: aRES × T7 gp5 co-folding
- [ ] T4: NTD-1 conservation / specificity-determining positions

Nothing in this repo contains invented sequences, structures or results. Uncertain values are marked
"VERIFY" in the data inventory.

## Repo layout

```
docs/                  the four documents above
data/
  sequences/           FASTA (record accession + date in SOURCE.md)
  structures/          PDB/mmCIF
  supplementary/       paper supplements (Tables S2/S3/S4 especially)
  external/            files uploaded or copied in from elsewhere
analysis/
  afold/               AlphaFold3 runs, one dir per run (T0.1, T1.1, ...)
  docking/             HADDOCK runs
  conservation/        alignments, HMMs, conservation scores, numbering maps
  scripts/             analysis code
results/               figures, tables, the construct list for the bench
refs/                  papers (gitignored — see below)
```

## Tooling

[`analysis/scripts/interface_contacts.py`](analysis/scripts/interface_contacts.py) — inter-chain
interface analysis for docking and co-folding output. Reads PDB **and** mmCIF (so both HADDOCK and
AlphaFold3 work), aggregates contact frequencies **across all models in a run** rather than trusting
one top-ranked pose, breaks the interface down by domain, and reports the P1 readout directly.
numpy only, no Biopython needed.

```bash
# Score a whole HADDOCK run against the hypothesis:
python3 analysis/scripts/interface_contacts.py analysis/docking/blind_run/*.pdb \
    --domain A:1-150:NTD1 --domain A:151-459:RES \
    --roi B:102-122:insertion1 \
    --tsv results/blind_contacts.tsv
```

The `--roi` line prints what fraction of the gp5-side interface falls in insertion1 — i.e. whether
the blind run independently supports P1. For AlphaFold models add `--min-plddt 70` to drop
low-confidence residues.

> Commit to scoring the blind run against P1 *before* looking at the output. A prediction registered
> in advance is worth considerably more than the same observation made afterwards.

## Note on `refs/`

PDFs are gitignored. The Osterman paper is a paywalled Elsevier article, and pushing it to GitHub
would redistribute it — if this repo is private and you want it tracked anyway, remove the
`refs/*.pdf` line from `.gitignore`.
