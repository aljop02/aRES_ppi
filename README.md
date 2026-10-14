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
| **[`docs/09_blind_zn_run_results.md`](docs/09_blind_zn_run_results.md)** | **The blind docking run with Zn modelled, scored against the pre-registered readout.** It did not converge — but the gp5 epitope came out blind, and one clean cluster puts acidic NTD-1 on the basic loop and yields a new two-mutant experiment |
| **[`docs/10_ortholog_zn_cofolds.md`](docs/10_ortholog_zn_cofolds.md)** | **Three aRES orthologs cofolded with Zn.** The Zn(Cys)₄ ribbon is confirmed with the metal present and conserved in all three. NTD-1 keeps its charge but is the fastest-diverging module; the charge *split* turns out to be UC5-only, and an apparent domain rearrangement is withdrawn as a PAE artifact |

`notes/vault/` is the Obsidian vault of prior work (notes + PyMOL/PAE figures), kept under version
control as the project record.

## The three things worth knowing up front

**1. The gp5 epitope is the 110–118 stretch, and recognition looks electrostatic.**
gp5 **S122P** and **D106G** both escape aRES (Fig. 4F), which places the epitope in the 100–125 loop.
The blind docking run then recovered that loop with no information supplied, and narrowed it:
**Arg111** and **Gly117** are contacted in **28/28** cluster representatives, while residues 104 and
105 are essentially absent (2/28 and 0/28) ([`docs/09` §2](docs/09_blind_zn_run_results.md)). The cleanest pose in the run is
charge-complementary — **acidic NTD-1** (net −16, zero lysines, and acidic in all three orthologs,
[`docs/10` §3](docs/10_ortholog_zn_cofolds.md)) against the loop's basic patch, with aRES
**Glu89/Glu93** salt-bridging gp5 **Arg111** ([`docs/09` §3](docs/09_blind_zn_run_results.md)).
>
> One caveat, and it is the load-bearing one: **that run did not converge**, and the loop is the most
> electrostatically attractive feature on gp5's surface, so a charge-dominated search should steer
> there regardless of biology. The epitope itself rests on the escape genetics; the charge model is a
> hypothesis, and point 3 is how it gets tested.

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
every gp5 and aRES point mutant is just site-directed mutagenesis. The set to run first:
**aRES E89A / E93A / D90A + gp5 R111A** — the first functional tests ever aimed at PF18870, and they
test the charge model rather than any one predicted contact
([`docs/09` §5](docs/09_blind_zn_run_results.md), [`docs/10` §6](docs/10_ortholog_zn_cofolds.md)).
Add **gp5 S122A**: S122P is a real escape mutation, so whether the serine hydroxyl or the proline's
backbone rigidity is what matters is a question the genetics already poses, and one primer pair
answers it. → [status §4](docs/05_status_after_vault_integration.md).

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
- [x] Blind docking run with Zn modelled (job 780832) — scored against the pre-registered readout in [`docs/09`](docs/09_blind_zn_run_results.md). **Did not converge**; Zn site verified correct; the pre-registered cofold contact not reproduced
- [x] Three orthologs cofolded with Zn → [`docs/10`](docs/10_ortholog_zn_cofolds.md). **Zn(Cys)₄ site confirmed with the metal present, same C8/C11/C31/C34 in 15/15 models**; NTD-1 acidic in all three but the fastest-diverging module; charge split is UC5-only
- [ ] **Record the two ortholog accessions** (marked VERIFY in [`SOURCE.md`](analysis/afold/cofold_zn_three_species/SOURCE.md)) and add both sequences to `data/sequences/` — `docs/10` §2 is unciteable without them
- [ ] **Fetch the unclustered HADDOCK score −154.0 structure** from the server — best pose in the run, absent from the summary archive, deleted after two weeks
- [ ] **AF3 seed replicates of the cofold** — as a general interface-reproducibility check: does *any* NTD-1 contact recur across seeds, and at what ipTM? Not to rescue a particular residue pair
- [ ] **aRES E89A / E93A + D90A + gp5 R111A** through the 4G/4H assays — E89A tests a conserved structural role, E93A a UC5-specific one, and they are not interchangeable ([`docs/10` §6](docs/10_ortholog_zn_cofolds.md))
- [ ] **10–20 PF18870+RES orthologs** through the same per-region identity analysis, with phylogenetic correction — pure compute, and the test that decides whether NTD-1's divergence is real ([`docs/10` §8](docs/10_ortholog_zn_cofolds.md))
- [ ] **Reconcile the UC5 frame** (model has Ile234/Ser293; paper says S234A/V293A) → then the UC5↔S_34AVR map
- [ ] **λ P / O / Gam co-expression with UC5 aRES** — tests the non-polymerase trigger
- [ ] **UC5 × GN04826 NTD-1 chimeras**, T7 readout (long lead time — start now). Add the two cofolded orthologs as donors: *Klebsiella* NTD-1 is the informative test of modularity, *E. chengduensis* the folding control ([`docs/10` §7](docs/10_ortholog_zn_cofolds.md))
- [ ] gp5 S122A through the 4G/4H assays — hydroxyl vs. proline rigidity, from the escape genetics alone
- [ ] T2 vs T4 gp43 alignment — zero cost, natural specificity pair
- [ ] Re-measure the three mutation panels with a consistent atom pair
- [ ] Supplementary Tables S2/S3/S4 acquired → conservation analysis

Nothing in this repo contains invented sequences, structures or results. Uncertain values are marked
"VERIFY" in the data inventory.

## Repo layout

```
docs/                  the documents listed above (01-10)
data/
  sequences/           FASTA (record accession + date in SOURCE.md)
  structures/          PDB/mmCIF
  supplementary/       paper supplements (Tables S2/S3/S4 especially)
  external/            files uploaded or copied in from elsewhere
analysis/
  afold/               AlphaFold3 runs, one dir per run (T0.1, T1.1, ...)
                       cofold_zn_three_species/  the three ortholog+Zn jobs (docs/10)
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
    --chains A B --domain B:1-120:NTD1 --domain B:126-447:RES \
    --roi A:100-125:loop --tsv results/blind_contacts.tsv
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

[`analysis/scripts/haddock_cluster_stats.py`](analysis/scripts/haddock_cluster_stats.py) —
recovers the full per-cluster statistics of a HADDOCK run from the server's Bokeh analysis plots.
The `*_summary.tgz` download gives only the top four structures per cluster, which **hides the two
numbers that decide whether a run converged**: the cluster populations, and whether the
best-scoring structure is clustered at all. Both are embedded in the plot HTML as base64 arrays.
Run it before HADDOCK deletes the results — afterwards the populations are unrecoverable.

```bash
python3 analysis/scripts/haddock_cluster_stats.py <run>/plots/fcc_hs.html \
    --tsv results/run_cluster_stats.tsv
```

In job 780832 this is what showed that the best pose in the run (HADDOCK score −154.0) was an
unclustered singleton, and that the two best-scoring clusters held 9 and 8 structures out of 207 —
neither visible from the summary archive.

[`analysis/scripts/af3_zn_site.py`](analysis/scripts/af3_zn_site.py) — the companion to
`zn_site_check.py` for AF3 jobs where **Zn was supplied as a ligand**, so the coordination shell can
be measured instead of inferred. Reports the shell, bond lengths against the 2.30 Å Zn–S reference,
tetrahedral RMSD, ligand and ion pLDDT, how often the same shell recurs across the five models, and
the job's confidence metrics. Only side-chain donor atoms define the geometry — a Cys CB sits ~2.8 Å
from the metal because its own SG is bonded to it, and counting it as a fifth ligand inflates the
tetrahedral RMSD from ~5° to ~28° on a perfect site.

```bash
python3 analysis/scripts/af3_zn_site.py <af3_output_dir> --cutoff 3.5 \
    --regions 1-48:Zn-ribbon 49-120:NTD1 --tsv results/zn_cofold_sites.tsv
```

[`analysis/scripts/af3_pae_blocks.py`](analysis/scripts/af3_pae_blocks.py) — block-averaged PAE from
the `full_data_*.json` files, with a verdict line per block pair. **Run this before making any
inter-domain claim from an AF3 model.** pLDDT is happily high for two domains whose relative
placement is a guess; PAE is the only number that licenses a statement about arrangement. In
`docs/10` §4 this is what retracted an apparent domain rearrangement: the N-block-to-RES PAE is
12.5–22.8 Å in all three orthologs, so the orientation was never determined and the cross-frame
RMSD was measuring AF3's ignorance.

```bash
python3 analysis/scripts/af3_pae_blocks.py <af3_output_dir> \
    --blocks 1-48:Zn-ribbon 49-120:NTD1 126-447:RES 448-459:tail
```

[`analysis/scripts/homolog_compare.py`](analysis/scripts/homolog_compare.py) — maps residue
numbering across orthologs and superposes them. Every residue number in this project is in the UC5
frame, and orthologs differ in length (*Klebsiella* aRES is 443 aa against UC5's 459), so quoting
"E93" for another ortholog means nothing until the frames are aligned — and eyeballing an alignment
is how frame errors get into documents. Global pairwise alignment (BLOSUM62, affine gaps), then
per-region identity with a two-proportion test, indel positions, invariant positions, and Kabsch CA
RMSD over chosen regions. Sequences are read from the structures themselves so the alignment and the
coordinates cannot drift apart.

```bash
python3 analysis/scripts/homolog_compare.py uc5.cif ortholog1.cif ortholog2.cif \
    --names UC5 E.chengduensis Klebsiella \
    --regions 1-48:Zn-ribbon 49-120:NTD1 126-447:RES \
    --map 8 11 31 34 89 90 93 --conserved 49-120:NTD1 --super 1-48:Zn-ribbon 49-120:NTD1
```

## Note on `refs/`

PDFs are gitignored. The Osterman paper is a paywalled Elsevier article, and pushing it to GitHub
would redistribute it — if this repo is private and you want it tracked anyway, remove the
`refs/*.pdf` line from `.gitignore`.
