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
| [`docs/02_hypothesis_NTD1_sensor.md`](docs/02_hypothesis_NTD1_sensor.md) | The hypothesis, the evidence for it, alternatives to exclude, numbered predictions P1–P9 |
| [`docs/03_computational_plan.md`](docs/03_computational_plan.md) | Tiered plan T0–T6, ordered by information per unit effort |
| [`docs/04_data_inventory.md`](docs/04_data_inventory.md) | Every sequence/structure/supplement needed, with accessions and acquisition status |
| [`docs/05_status_after_vault_integration.md`](docs/05_status_after_vault_integration.md) | Reconciles `01`–`04` with the prior work in `notes/vault/`: real domain boundaries, the cofold result, and what the earlier docs got wrong |
| **[`docs/06_phage_range_matrix.md`](docs/06_phage_range_matrix.md)** | **Current state.** Fig. 1B decoded, the chimera pair fixed, and the λ problem that breaks the polymerase model. **Read this one first if you only read one.** |
| [`docs/07_docking_run_review.md`](docs/07_docking_run_review.md) | The Zn(Cys)₄ ribbon found in aRES, review of the blind-run parameters, and where blind vs restrained docking each help |
| **[`docs/08_embl_hamburg_structural_plan.md`](docs/08_embl_hamburg_structural_plan.md)** | **The rotation plan for EMBL Hamburg.** Two tracks (first PF18870 structure / the interaction), the gate structure, week 1 in full, what was dropped and why, which negative results mean nothing, and the facility claims to verify on arrival |
| **[`docs/09_haddock_blind_run_results.md`](docs/09_haddock_blind_run_results.md)** | **All three HADDOCK runs.** The two blind runs scored against P1 as pre-registered — neither converged, so P1 is neither supported nor falsified, and the score is an electrostatic artifact traceable to the unmodelled Zn site. Plus §10: the restrained run `2_ares_T7` failed outright with an empty restraint table, and how to stop that happening again |

`notes/vault/` is the Obsidian vault of prior work (notes + PyMOL/PAE figures), kept under version
control as the project record.

## The three things worth knowing up front

**1. The two escape mutations are one mechanism, and the cofold already supports it.**
gp5 **S122P** and **D106G** both escape aRES. In the existing cofold, gp5 **Ser122 OG** hydrogen-bonds
to aRES **Asp104** (3.06 Å) — and Asp104 sits inside the NTD-1 block. Meanwhile **Asp106** caps the
**Gly121** backbone amide *within* gp5 (2.71 Å). So: **D106 is the scaffold that presents the loop,
S122 is the contact.** That explains why two chemically unrelated mutations 16 residues apart give the
same phenotype. → [status §2](docs/05_status_after_vault_integration.md).

**2. T7 DNA polymerase is a heterodimer, and the activating species was the holoenzyme.**
gp5 (704 aa) binds host *E. coli* **thioredoxin** as its processivity factor. The enzyme used in the
paper's in vitro activation experiment (NEB M0274S) is the **gp5·TrxA complex**, not gp5 alone.
Docking gp5 by itself may be docking the wrong species.
→ [hypothesis §4](docs/02_hypothesis_NTD1_sensor.md).

**2b. λ-vir has no DNA polymerase — yet two of the three aRES proteins defend against it.**
λ replicates with host Pol III; its only replication proteins are O and P. So **at least one aRES
trigger is not a polymerase**, and the paper's framing holds for T7 but is not the general rule.
Possible resolution: NTD-1 may sense phage proteins that *hijack host replication factors* —
T7 gp5 captures thioredoxin, λ P captures DnaB. Testable in one co-expression assay.
→ [matrix §2](docs/06_phage_range_matrix.md).

**3. The decisive experiments are cheap and need no mutant phages.**
The paper's co-expression toxicity assay (Fig. 4G) and pull-down (Fig. 4H) are both plasmid-based, so
every gp5 and aRES point mutant is just site-directed mutagenesis. Two to run first:
**gp5 S122A** (removes the hydroxyl donor without proline's backbone rigidity → separates direct
contact from loop conformation) and **aRES D104A** (the reciprocal test, and the first functional
test ever aimed at PF18870). → [status §4](docs/05_status_after_vault_integration.md).

Longer lead time but the cleanest test of the central hypothesis: swap NTD-1 between **UC5 and
GN04826**, which have near-identical phage ranges but differ by ~5 logs at T7 alone — so the other
eight phages act as built-in folding controls in the same assay (P2,
[matrix §1](docs/06_phage_range_matrix.md)). **Start that cloning in parallel, not after.**

## Status

- [x] Paper read and summarized
- [x] Hypothesis formalized with falsifiable predictions
- [x] Computational plan drafted
- [x] Interface-analysis tooling written and tested
- [x] Prior work integrated (`notes/vault/`): domain boundaries, aRES×gp5 cofold, MD plan
- [x] Fig. 1B decoded — GN04826 confirmed **not** to defend against T7; chimera pair fixed
- [x] **Blind HADDOCK runs 779823 + 780456 scored against P1** (`docs/09`) — neither converged
      (~94% unclustered); P1 neither supported nor falsified; scoring dominated by an electrostatic
      artifact from the unmodelled Zn(Cys)₄ site
- [ ] **Download the two top-scoring models before the server expires them** — both are unclustered
      singletons, so they are *not* in the cluster-summary archives (`docs/09` §8)
- [ ] **Resubmit `2_ares_T7` with the explicit `ambig.tbl`** — the first attempt (job 780457) came
      back with an empty restraint table, so the partners never touched: 400 structures with the two
      chains ~270 Å apart (`docs/09` §10). Restraints are now generated by
      `analysis/scripts/make_ambig_tbl.py`; run a 100-structure test job and confirm Eair ≠ 0 and
      BSA > 0 before committing the full sampling
- [ ] **Re-dock with Zn²⁺ modelled** into the Zn(Cys)₄ ribbon, or exclude aRES 1–48
- [ ] **Dock gp5·TrxA** (T1.3) — 1T7P with the DNA stripped; gp5-alone docking cannot read the TBD
      diagnostic in `docs/02` §4
- [ ] **Reconcile the UC5 frame** (model has Ile234/Ser293; paper says S234A/V293A) → then the UC5↔S_34AVR map
- [ ] **λ P / O / Gam co-expression with UC5 aRES** — tests the non-polymerase trigger
- [ ] **UC5 × GN04826 NTD-1 chimeras**, T7 readout (long lead time — start now)
- [ ] gp5 S122A + aRES D104A through the 4G/4H assays
- [ ] T2 vs T4 gp43 alignment — zero cost, natural specificity pair
- [ ] Re-measure the three mutation panels with a consistent atom pair
- [ ] Supplementary Tables S2/S3/S4 acquired → conservation analysis

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
# In the finished HADDOCK runs, chain A is aRES and chain B is gp5 (see docs/09).
python3 analysis/scripts/interface_contacts.py analysis/docking/blind_run/*.pdb \
    --chains A B \
    --domain A:1-48:Zn_ribbon --domain A:49-120:NTD1 --domain A:126-447:RES \
    --roi B:102-122:gp5_insertion1 --tsv results/blind_contacts.tsv
```

The domain breakdown says whether the aRES side of the interface is NTD-1 or RES; the `--roi` line
says what fraction of the gp5 side falls in the epitope loop. For AlphaFold models add
`--min-plddt 70` to drop low-confidence residues.

> Commit to scoring the blind run against P1 *before* looking at the output. A prediction registered
> in advance is worth considerably more than the same observation made afterwards.

[`analysis/scripts/zn_site_check.py`](analysis/scripts/zn_site_check.py) — detects structural
Zn(Cys)₄ sites in models that contain no metal, by testing the SG···SG geometry against the ideal
tetrahedral edge (3.76 Å). Also reports mean pLDDT per region, which is how you decide whether a
loop or tail is worth docking.

```bash
python3 analysis/scripts/zn_site_check.py ares_model.pdb \
    --regions 1-48:Zn-ribbon 49-120:NTD1 448-459:tail
```

## Note on `refs/`

PDFs are gitignored. The Osterman paper is a paywalled Elsevier article, and pushing it to GitHub
would redistribute it — if this repo is private and you want it tracked anyway, remove the
`refs/*.pdf` line from `.gitignore`.
