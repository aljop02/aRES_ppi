# Does NTD-1 diverge fastest? Eighteen orthologs, phylogenetically corrected

Panel and provenance: [`data/sequences/SOURCE.md`](../data/sequences/SOURCE.md).
Machine-readable output: [`results/panel18_region_identity.tsv`](../results/panel18_region_identity.tsv),
[`results/panel18_tests.tsv`](../results/panel18_tests.tsv),
[`results/panel18_charge.tsv`](../results/panel18_charge.tsv),
[`results/pool190_tests.tsv`](../results/pool190_tests.tsv),
[`results/pool190_zn_ligand_map.tsv`](../results/pool190_zn_ligand_map.tsv),
[`results/pool190_pfam_bounds_in_uc5_frame.tsv`](../results/pool190_pfam_bounds_in_uc5_frame.tsv).

**Verdict against the three outcomes registered in [`docs/10`](10_ortholog_zn_cofolds.md) §8: the
third one matched. The Zn ribbon diverges as fast as NTD-1 — the two are statistically
indistinguishable (PGLS p = 0.36) — so the bimodular reading of the N-block in
[`docs/07`](07_docking_run_review.md) §1 is wrong, and the N-block is one evolutionary unit, not
two.** What *is* significant is the N-block against the catalytic domain: NTD-1 is 6.2 points less
conserved than the RES side (PGLS p = 0.0022, Holm-adjusted 0.0044) and 10.6 points less than the
PF08808 core (0.0004, adjusted 0.0008), and the ribbon is below the RES core too (0.0020, adjusted
0.0040). Every p-value in this document is phylogenetically corrected unless it says otherwise, and
the Holm adjustment for the several region pairs is quoted alongside. So `docs/10` §2's
direction of travel survives, but its subject changes: **the fast-diverging module is the whole
N-block, not NTD-1 within it.**

The registered outcome that did *not* match is the one the project wanted. Outcome 1 required
NTD-1 "significantly below **both** flanking modules"; it is below one of them. `docs/10` §2 can
now be cited for "the N-block diverges faster than the NADase it activates" and **cannot** be cited
for "NTD-1 is the fastest-diverging module", which is what it currently says.

**Separately, and more strongly: the charge result from `docs/10` §3 is the one that passes
decisively.** NTD-1 is net-acidic in **18/18** panel orthologs and in **187/190** of the wider set,
never net-positive in any of the 190, at −16.5 charges per 100 residues against −1.2 for the RES
core across those 190; on the 18-member panel that difference is PGLS p = 1×10⁻⁹, with the caveat
about such exponents in §3. That is a far sharper family-level signal than identity, exactly as
`docs/10` §8 predicted it would be. The specific counts quoted there do **not** generalise: "zero
lysines" holds in 64/190 and "exactly eleven aspartates" in neither panel nor pool.

No trigger protein was in any of this. It is sequence analysis of a predicted gene family and says
nothing directly about phage-protein recognition.

---

## 0. What was run

A panel of **18 PF18870+PF08808 proteins**, 25.3–41.0% identical to UC5, one per genus, spanning
18 genera and 4 phyla. Selection filters were fixed before any per-region identity was computed and
are listed in [`SOURCE.md`](../data/sequences/SOURCE.md); they condition on overall divergence and
on alignment coverage, never on the identity of one region against another.

Two things about the family are worth knowing before the numbers:

- **213 of the 233 PF18870 proteins in UniProtKB (91%) carry a RES domain.** The fusion is the
  normal state of this domain, not an *E. coli* peculiarity.
- **UC5's closest available relative is 41% identical.** There is nothing in UniProtKB at the 70%
  level where `docs/10`'s *E. chengduensis* sat — indeed no *Escherichia* or *Klebsiella* protein
  carries PF18870 at all. The panel is unavoidably all-distant, which costs statistical power and
  puts the alignments near the twilight zone. §6 reports what was done about that.
- **All 18, and all 213, are UniProtKB *unreviewed* entries with protein existence "Predicted".**
  None has experimental evidence of expression or of anti-phage function.

`docs/10`'s two orthologs are **not** in this panel — the filters were fixed on the annotated set
before they could be identified — but they are **no longer VERIFY.** Both sequences turned out to be
recoverable from the AF3 job archives (an AF3 model carries its own sequence, which is why
`homolog_compare.py` reads structures at all), and the CRC64 → UniParc route then identified all
three of the repo's unaccessioned sequences, verified position by position rather than by checksum
alone ([`SOURCE.md`](../data/sequences/SOURCE.md) has the method and the dead-accession list):

| | UniParc | RefSeq | organism |
|---|---|---|---|
| UC5 aRES | UPI00038F43E8 | WP_021575732 | *E. coli* — the same 459-aa sequence in ~150 assemblies |
| `docs/10` Kleb_aRES | UPI0006936E4D | WP_052450931 | *Klebsiella pneumoniae* — label confirmed |
| `docs/10` aRES_Echengduensis | UPI0003BED449 | WP_023337974 | ***Enterobacter*** (*E. cloacae* complex) |

Three things worth carrying forward. **Quote RefSeq, not UniProt**: every UniProtKB entry for all
three is deleted, including six for UC5. **`docs/10`'s "E. chengduensis" is an *Enterobacter***, so
the label is right if "E." meant *Enterobacter* and wrong if it meant *Escherichia*; ⚠️ the species
itself is not recoverable, because the identical sequence occurs across several *E. cloacae*-complex
members. And the deletions are the reason UC5 has no PF18870-annotated relative in UniProtKB: the
redundant *E. coli* entries were removed, not the protein — see §6.

### Per-region identity to UC5

Global pairwise alignment (BLOSUM62, affine gaps), regions in the UC5 frame:

| region | panel mean ± sd | range | median columns | mean coverage |
|---|---|---|---|---|
| Zn ribbon 1–48 | 27.9 ± 5.9 | 17.0–39.6 | 47 | 98% |
| **NTD-1 49–120** | **25.6 ± 4.6** | 19.4–39.4 | 65 | 89% |
| N-block 1–120 | 26.5 ± 4.1 | 21.1–39.5 | 112 | 93% |
| RES side 126–447 | 31.7 ± 4.4 | 25.1–41.7 | 243 | 80% |
| **PF08808 core 199–363** | **36.2 ± 8.3** | 23.2–55.5 | 148 | 89% |

"Median columns" is what each ortholog's own number rests on. A 47-column estimate and a
243-column estimate are not comparable in precision, which is the first reason the test below is
not a test on positions.

---

## 1. 🔬 The test, and what the corrections did to it

`docs/10` §2's p = 0.037 came from a two-proportion z-test on alignment *positions* in **one**
ortholog. That treats 71 positions as 71 independent draws and one ortholog as a sample. Both
assumptions are replaced here:

- the unit of analysis is the **ortholog**, not the position — one paired difference
  *d* = identity(NTD-1) − identity(other) per ortholog, so the position count leaves the test and
  survives only as the precision estimate above;
- **Holm–Bonferroni** across the region pairs compared, because comparing three regions is several
  tests and `docs/10` §2 quoted the best of them;
- **PGLS** against a Brownian covariance built from the panel's own all-against-all identities, so
  that orthologs sharing ancestry do not count as independent observations.

The registered comparison, exactly as `docs/10` §8 framed it — NTD-1 against both flanking
modules, two tests:

| | mean *d* | orthologs with *d* < 0 | Wilcoxon (exact) | Holm | PGLS | Holm |
|---|---|---|---|---|---|---|
| **NTD-1 − RES side** | **−6.16** | **17/18** | 7.6×10⁻⁵ | 1.5×10⁻⁴ | **0.0022** | **0.0044** |
| NTD-1 − Zn ribbon | −2.30 | 12/17 | 0.159 | 0.159 | 0.364 | 0.364 |

The ribbon row counts 17, not 18, because *Mycobacterium senriense* `A0ABM7SVF0` is 25.0% identical
in both regions — an exact tie, which Wilcoxon's and the sign test's own rules drop.

Wilcoxon signed-rank and sign tests are exact (enumerated, not normal-approximated); the paired
*t*-test agrees with both throughout and is in
[`results/panel18_tests.tsv`](../results/panel18_tests.tsv). The PGLS column is the
Brownian-covariance p-distance model; the Poisson-corrected distance gives 0.0098 and 0.527, i.e.
the same verdict with less power. **Nothing here is significant for the ribbon comparison under any
test, corrected or not.**

### Supplementary: the same question with a coverage-matched comparator

"RES side 126–447" is a poor comparator for a reason worth recording. The UC5 segment **364–447 is
aligned in only 53% of the panel on average (range 14–100%)**, against 89% for 199–363 and 98% for
the ribbon — most orthologs simply do not have it. The PF08808 core is the matched comparator:

| | mean *d* | *d* < 0 | PGLS | PGLS, indels counted |
|---|---|---|---|---|
| NTD-1 − PF08808 core | −10.60 | 17/18 | **0.0004** | **0.0012** |
| Zn ribbon − PF08808 core | −8.30 | 16/18 | **0.0020** | — |
| N-block 1–120 − PF08808 core | −9.67 | 17/18 | **0.0003** | — |
| NTD-1 − Zn ribbon | −2.30 | 12/17 | 0.364 | 0.054 |

Against the matched comparator every N-block-vs-catalytic-domain comparison is significant and the
NTD-1-vs-ribbon comparison still is not. **This is the cleanest statement the data support: both
N-terminal modules diverge faster than the PF08808 core, at indistinguishable rates from each
other.**

### Where the two identity definitions disagree, and why

`id_aligned` divides by aligned columns, so a deletion inside a region costs nothing; `id_allpos`
divides by every UC5 position, so it costs. Against the unmatched `RES 126–447` the choice flips
the headline result — NTD-1 − RES goes from p = 0.0022 to p = 0.15 — and the flip is entirely the
missing 364–447: counting absent positions as mismatches punishes the region that is absent.
Against the coverage-matched core the two definitions agree (p = 0.0004 and 0.0012). **The headline
rests on the matched comparator, where the metric choice does not matter.** Quoting the
`126–447` number without that caveat would be quoting an artifact of region definition.

---

## 2. What the phylogenetic correction actually bought

Re-run on the full **190-sequence** qualifying set — the panel before genus de-duplication, which
contains sequences up to 100% identical to each other:

| | mean *d* | *d* < 0 | Wilcoxon, uncorrected | PGLS | PGLS, Poisson |
|---|---|---|---|---|---|
| NTD-1 − RES side | −5.91 | 153/190 | < 10⁻⁴ | **< 10⁻⁴** | **0.0037** |
| NTD-1 − Zn ribbon | −3.05 | 122/189 | **< 10⁻⁴** | **0.095** | **0.304** |

**The second row is the whole argument for doing this properly.** At n = 190 the uncorrected test
reports p < 0.0001 for the ribbon comparison and the phylogenetic correction puts the same data at
p = 0.09–0.30. Adding 172 sequences made the naive p-value smaller and changed nothing about the
evidence, because the added sequences were mostly near-relatives of ones already there. Had
`docs/10` §2 simply been scaled up without correction, it would have reported a highly significant
result for precisely the comparison that does not hold.

The RES comparison, by contrast, is significant with and without correction, at n = 18 and at
n = 190, under both distance models. That is the part to trust.

---

## 3. 🔬 The charge result is the strong one

`docs/10` §3 proposed this as "a sharper test than identity". It is.

Net charge (Lys+Arg − Asp−Glu) per 100 residues, regions mapped onto each ortholog **through the
alignment** rather than sliced out of its own numbering — an ortholog with a 14-residue insertion
at UC5 ref~66 does not have NTD-1 at 49–120 in its own frame:

| region | panel (n=18) | pool (n=190) | net < 0, pool |
|---|---|---|---|
| Zn ribbon | −6.0 ± 7.3 | −5.9 ± 7.9 | 132/190 (69%) |
| **NTD-1** | **−17.7 ± 3.3** | **−16.5 ± 4.8** | **187/190 (98%)** |
| N-block | −11.9 ± 4.4 | −11.9 ± 4.6 | 185/190 (97%) |
| RES side | −0.7 ± 3.1 | −0.9 ± 3.2 | 126/190 (66%) |
| PF08808 core | −0.2 ± 3.0 | −1.2 ± 2.9 | 120/190 (63%) |

**NTD-1 is never net-positive in 190 orthologs.** The three exceptions to "net < 0" are exactly
neutral, not positive. The other regions sit near 63–69%, i.e. close to what a region of no
particular charge character does. Paired, corrected, on the panel:

| | mean *d* | *d* < 0 | PGLS | Holm |
|---|---|---|---|---|
| NTD-1 − Zn ribbon | −11.62 | 16/18 | 2.0×10⁻⁴ | 2.0×10⁻⁴ |
| NTD-1 − RES side | −16.95 | **18/18** | 1.2×10⁻⁸ | 2.4×10⁻⁸ |
| NTD-1 − PF08808 core | −17.47 | **18/18** | 1.0×10⁻⁹ | 3.1×10⁻⁹ |

Two honest limits on those exponents. **The exact Wilcoxon test cannot resolve below
2/2¹⁸ = 7.6×10⁻⁶** on 18 paired observations, and the 18/18 rows sit exactly at that floor — the
10⁻⁸ and 10⁻⁹ values are the parametric PGLS *t*-test extrapolating past where the
distribution-free test runs out, so read them as "as strong as 18 orthologs can show", not as
literal probabilities. And charge composition is a far smoother trait than identity; a module could
be acidic because of what it does chemically rather than because of selection on the module as a
unit.

**Unlike identity, this separates NTD-1 from the ribbon** (p = 2×10⁻⁴). So the two N-terminal
modules evolve at the same *rate* while remaining chemically distinct — which is what two
co-evolving halves of one domain with different jobs should look like, and is consistent with
`docs/10` §1's PAE finding that they are one structural module.

### What does not survive

- **"Zero lysines in NTD-1" is not a family feature.** 64/190 (34%) have none; the panel range is
  0–4 lysines and the pool's is 0–12. It is a real enrichment — 2/190 for the RES side, 8/190 for
  the core — but it is not an invariant, and `docs/08` §3's "no lysine anywhere in 49–120" is a
  UC5 fact being used as a design constraint, which is fine, as long as it is not read as a family
  property. It bears directly on `docs/08`'s P2 XL-MS experiment, which needs **gp5** lysines, not
  aRES ones, so that design is unaffected.
- **"Exactly eleven aspartates" is a coincidence of three sequences.** Panel range 4–15, pool 3–17.
- **The acidity is in density, not in counts.** Per 100 residues NTD-1 spans −26.0 to 0.0 across
  190 orthologs; that is the number to quote.

---

## 4. The Zn site: the positions hold family-wide, "Cys₄" does not

`docs/10` §1 concluded that `docs/07`'s "never mutate C8/C11/C31/C34" "now holds family-wide".
Mapping those four UC5 positions onto all 190 orthologs
([`results/pool190_zn_ligand_map.tsv`](../results/pool190_zn_ligand_map.tsv)):

| UC5 position | Cys | His | other / gap |
|---|---|---|---|
| C8 | 185 (97%) | — | 1 Gly, 4 gap |
| **C11** | **145 (76%)** | **40 (21%)** | 1 Ala, 4 gap |
| C31 | 189 (99%) | — | 1 Tyr |
| C34 | 189 (99%) | — | 1 Ala |

| shell | n | % |
|---|---|---|
| **Cys₄ (CCCC)** | **144** | **76%** |
| **Cys₃His (CHCC)** | **39** | **21%** |
| anything else | 7 | 4% |

**All four positions are Cys or His in 183/190 (96%), but the shell is Cys₄ in only 76% — in a
fifth of the family position 11 is histidine.** The mutagenesis advice is unchanged and if anything
reinforced: these are the most positionally conserved residues in the protein, more conserved than
anything in the PF08808 core. But the site should be described as a **Zn ribbon with a Cys₃His
variant**, not as Cys₄, and two consequences follow:

- a construct designed around "four cysteines compete with six histidines for Zn(II)"
  (`docs/08` §5, item 7's reason for C-terminal tagging) is sound for UC5, and the His variant does not
  change that, since the variant histidine is a ligand rather than a competitor;
- the homolog-expansion panel in `docs/08` §5, item 7 should **deliberately include one Cys₃His ortholog**.
  If the ribbon's job is metal binding, it should fold; if the project's reading of a strict Cys₄
  ribbon is load-bearing anywhere, that is the cheap test of it.

This also sharpens §1's result. The ribbon's sequence identity is low (27.9%, no better than
NTD-1) while its ligand positions are near-invariant — low identity across a module is not evidence
of low constraint on the residues that matter in it.

---

## 5. ⚠️ PF18870's own boundaries contradict the repo's `49–136`

This is a side finding of the mapping work, and it matters for `docs/08`.

The repo records **InterPro PF18870 = 49–136** for UC5 in
[`docs/02`](02_hypothesis_NTD1_sensor.md), [`docs/03`](03_computational_plan.md) §T1,
[`docs/05`](05_status_after_vault_integration.md) §1.1, [`docs/07`](07_docking_run_review.md) §1 and
[`docs/08`](08_embl_hamburg_structural_plan.md) §3, and `docs/08` §7/D1 builds an argument on it:
"PF18870's annotated start at 49 is exactly where the Zn ribbon ends, so Pfam is blind to the
ribbon rather than reporting a folding boundary", with `hmmalign` against PF18870 said to be unable
to align 1–48 for that reason (`docs/08` §5, item 2).

Mapping each ortholog's **own InterPro-annotated** Pfam bounds back into the UC5 frame through the
alignment ([`results/pool190_pfam_bounds_in_uc5_frame.tsv`](../results/pool190_pfam_bounds_in_uc5_frame.tsv)):

| | median in UC5 frame | IQR | range | repo's figure |
|---|---|---|---|---|
| PF18870 | **25–147** | 25–27 / 144–149 | 23–41 / 125–162 | 49–136 |
| PF08808 | **201–365** | 181–207 / 357–394 | 179–212 / 352–444 | 199–363 ✅ |

- **PF18870 starts before UC5 49 in 190/190 orthologs**, and at or before residue 34 in 189/190 —
  i.e. it covers the ribbon's *second* CxxC (C31/C34) but not its first (C8/C11).
- **It ends after UC5 120 in 190/190**, median 147.
- **PF08808 maps to 201–365, matching the repo's 199–363 to within 2 residues** — so the method
  reproduces one of the two numbers it is being checked against, which is why the PF18870
  discrepancy is worth taking seriously rather than dismissing as alignment noise.

So the family-level claim is wrong: **Pfam is not blind to the Zn ribbon.** The domain boundary
falls *between the two CxxC motifs*, which is a more interesting place for it to fall and is
consistent with everything in §1 and §3 — the Pfam domain, the AF3 PAE module
(`docs/10` §1: ribbon↔NTD-1 PAE 5.0–6.4 Å) and the evolutionary unit found here all say the
N-block is one piece.

⚠️ **VERIFY, and this one is cheap and load-bearing.** The UC5-specific number could not be
re-derived here: UC5 aRES is not in UniProtKB, no HMMER is installed locally, and the InterProScan
REST API requires an email address which was not sent. Two readings remain open — the repo's
`49–136` came from an InterProScan run that genuinely matched UC5 only from 49, or it is an
envelope-versus-aligned-coordinate mix-up or an older Pfam release. **Run `hmmscan` with the
PF18870 HMM against `data/sequences/ares_uc5.fasta` and settle it before ordering the boundary
series.** `docs/08` §7/D1 already flags that "1–120 truncates the Pfam domain by 16
residues"; if the orthologs are right, 1–120 truncates it by **27**, and the `49–136`
"no-ribbon" construct is not a no-ribbon construct — it would start 24 residues inside the domain
and cut C31/C34 off from the rest of it. That changes what a solubility failure of `49–136` would
mean.

I have not edited `docs/02`, `03`, `05`, `07` or `08`; the number is theirs to correct once the
HMM call is in.

---

## 6. Robustness: what was checked, and what it changed

| check | NTD-1 − RES | NTD-1 − ribbon |
|---|---|---|
| panel of 18, PGLS | **0.0022** | 0.364 |
| full 190, PGLS | **< 10⁻⁴** | 0.095 |
| Poisson-corrected distances | **0.0098** | 0.527 |
| collapse sequences ≥40% identical (14 clusters) | **0.0006** | 0.147 |
| collapse ≥50% (17 clusters) | **0.0002** | 0.193 |
| leave one ortholog out, worst of 18 | **0.0002** | 0.274 |
| leave one phylum out, worst of 4 | **0.039** (drops 10 of 18) | 0.636 |
| indels counted as mismatches, matched comparator | **0.0012** | 0.054 |

The RES comparison survives every one of these. The ribbon comparison fails all of them. The
weakest RES result is dropping Pseudomonadota, which removes 10 of 18 orthologs and leaves
p = 0.039 — a real weakening, honestly the limit of an 18-member panel, and the reason the 190-set
re-run matters: there the same jackknife leaves p < 10⁻⁴ with 105 sequences remaining.

**What was not done.** There is no gene tree: the covariance comes from a UPGMA tree built from the
same pairwise identities being tested, which is circular to a degree I cannot quantify, and the
trait (identity *to UC5*) is not a per-species character in the strict comparative sense. A proper
treatment would build a maximum-likelihood tree from a multiple alignment and fit per-branch rates
per region. That is the next increment if this ever needs to be publication-grade; it is unlikely to
rescue the ribbon comparison, which fails at p = 0.36 with the correction *and* p = 0.16 without it.

**The panel is biased by how UniProtKB handles redundancy, and this is the limitation I would
flag hardest.** UC5's own sequence occurs in ~150 *E. coli* assemblies and has six UniProtKB
accessions, *all deleted* as redundant — which is why §0's search found no PF18870 protein in
*Escherichia* at all. The same process must thin every well-sampled genus. So the panel is drawn
from a database that preferentially retains proteins from under-sampled lineages, which plausibly
shifts it towards the distant end of the family and is part of why nothing above 41% was available.
It does not bias NTD-1 against RES *within* a protein — the test is paired, and both regions come
from the same sequence — but it does mean "the family" here is the family as UniProtKB chose to
represent it. A panel built from RefSeq WP\_ records, which are not thinned the same way, would be
the fix, and it would probably recover the 50–70% band that `docs/10`'s two orthologs occupy.

**Alignment reliability is the floor under all of this.** At 25–41% identity, pairwise global
alignment is near the twilight zone, and a region that is harder to align will show spuriously low
identity. Three things argue the effect is not driving the result: the ribbon, which takes the
fewest indels and aligns best (98% coverage), is *not* the most conserved region; the result holds
under the indel-inclusive metric against the matched comparator; and it holds in the
highest-identity orthologs as well as the lowest (no trend in the per-ortholog differences,
[`results/panel18_paired_diffs.tsv`](../results/panel18_paired_diffs.tsv)). A profile-based
alignment would still be the better instrument.

### Tooling

[`analysis/scripts/homolog_compare.py`](../analysis/scripts/homolog_compare.py) was extended rather
than replaced: it now takes FASTA as well as structures (multi-record, so a panel is one argument),
writes per-region identity counts for the statistics, and fills its alignment matrices a row at a
time instead of a cell at a time — **~100× faster, and verified to return byte-identical alignments
to the committed version on 31 test cases** including self-, reversed-, single-residue and
unknown-residue inputs. That mattered: the 190-sequence tree needs 17 955 pairwise alignments.

[`analysis/scripts/ortholog_panel_stats.py`](../analysis/scripts/ortholog_panel_stats.py) is new and
holds the statistics. No scipy in this environment, so every distribution is hand-rolled and
**checked against 38 analytic or exhaustively-enumerated targets** — exact Wilcoxon against brute
force over all 2ⁿ sign assignments including the tied case, the sign test against binomial
coefficients, Student-*t* against standard critical values and the Cauchy closed form, UPGMA against
a textbook five-taxon example, PGLS against the one-sample *t*-test when V = I, and Holm against
hand-computed adjustments.

[`analysis/scripts/charge_composition.py`](../analysis/scripts/charge_composition.py) gained panel
handling and `--ref`, which maps region bounds through an alignment instead of slicing the
reference's coordinates out of a different protein. It also now **refuses** a multi-record FASTA on
the single-sequence path, which previously concatenated every record into one sequence and reported
composition for a chimera that does not exist.

```bash
python3 analysis/scripts/homolog_compare.py data/sequences/ares_uc5.fasta \
    data/sequences/pf18870_res_panel18.fasta --names UC5 \
    --regions 1-48:Zn-ribbon 49-120:NTD1 126-447:RES \
    --identity-tsv results/panel18_region_identity.tsv

python3 analysis/scripts/ortholog_panel_stats.py results/panel18_region_identity.tsv \
    --panel data/sequences/pf18870_res_panel18.fasta \
    --focal NTD1 --against Zn-ribbon RES \
    --groups results/panel18_groups.tsv --out-prefix results/panel18
```

---

## 7. What this does and does not license

**Licensed:**
- The **N-block (1–120) diverges faster than the PF08808 core**, by ~10 identity points, in 17/18
  orthologs, phylogenetically corrected and robust to every check in §6 (§1).
- **NTD-1 and the Zn ribbon diverge at indistinguishable rates** (§1) — stated as a failure to
  distinguish on 18 orthologs, not as proof of equality.
- **NTD-1 is net-acidic family-wide**: 187/190 negative, never positive, and significantly more
  acidic than both flanking regions (§3).
- The **four Zn ligand positions are the most conserved positions in the protein**, with a Cys₃His
  variant in 21% of the family (§4).
- **PF18870 covers the ribbon's second CxxC and extends past UC5 120** in 190/190 orthologs (§5).

**Not licensed:**
- **"NTD-1 is the fastest-diverging module."** It is not distinguishable from the ribbon. `docs/10`
  §2's wording and the README's summary of it both need correcting.
- The **bimodular N-block** of `docs/07` §1, as an *evolutionary* claim. It may still be right
  structurally — two modules can be one evolutionary unit — but §1 is evidence against treating
  them as separately evolving, and §5 removes the Pfam-boundary argument for the split.
- **"Zero lysines" or "eleven aspartates" as family features** (§3).
- **The UC5 PF18870 boundary**, in either direction, until an HMM is run (§5).
- Anything about **phage-protein recognition**. No trigger was in this analysis, and a divergence
  rate is not a binding site. The modular-sensor hypothesis gains a *comparative* leg here, for the
  N-block rather than for NTD-1; it is not tested.

### What this means for the hypothesis

Per `docs/10` §8, outcome 3 means the N-block is one evolutionary unit. The modular-sensor
hypothesis is **not refuted and is mildly supported** — the recognition block does diverge faster
than the catalytic domain it activates, which is what a specificity module should do. But the unit
of specificity is now the whole N-block, and that has a concrete consequence: **the chimera
experiment should swap 1–120 (or 1–147, pending §5), not 49–120.** Swapping 49–120 would cut
through the middle of the Pfam domain, through the AF3 PAE module, and through the evolutionary
unit, all three of which agree on where the boundary is not. `docs/08` §7/D1 already treats
"1–120 soluble while 49–120 is not" as the informative comparison; this says the same thing from
sequence.

---

## 8. Next steps, cheapest first

1. **`hmmscan` the PF18870 HMM against `data/sequences/ares_uc5.fasta`** and settle §5. One command
   once HMMER is available, and it gates the `docs/08` boundary series. Until then the `49–136`
   figure in `docs/02`, `03`, `05`, `07` and `08` is VERIFY.
2. **Correct `docs/10` §2 and the README** from "NTD-1 is the fastest-diverging module" to "the
   N-block diverges faster than the RES core" — the README is updated here; `docs/10` §2 is left
   for whoever owns it, since another pass was editing that file while this ran.
3. **Add one Cys₃His ortholog to the `docs/08` §5 homolog-expansion order** (§4). Zero marginal
   cost in a gene order that already includes 3–4 orthologs.
4. **Swap 1–120, not 49–120, in the chimera panel** (§7).
5. **Rebuild the panel from RefSeq `WP_` records instead of UniProtKB** (§6). UniProtKB's
   redundancy removal is why nothing above 41% identity was available, and `docs/10`'s own two
   orthologs sit at 44.5% and 71.9% — the informative band is missing from the panel for a database
   reason, not a biological one. This is the single change most likely to sharpen §1, because the
   NTD-1-vs-ribbon comparison fails on power rather than on direction.
6. **A profile alignment and an ML gene tree**, if §1 ever needs to be publication-grade (§6). The
   honest version of this analysis is a per-region rate model on a real tree; what is here is the
   correction that the available tooling supports, and it is enough to decide between the three
   registered outcomes.
7. **The conserved positions *within* NTD-1 are the remaining mutagenesis question.** §1 says the
   module diverges; §4 shows that a diverging module can still hold near-invariant positions. The
   per-position invariance profile across 190 orthologs is one `--conserved` run away and is the
   thing that would actually nominate residues.

### Registered in advance, for item 7

- **NTD-1 contains a set of near-invariant positions comparable in density to the ribbon's ligand
  set** → those positions, not the variable ones, are the mutagenesis targets, and the
  "specificity module" reading has to accommodate a conserved core.
- **NTD-1's invariant positions are sparse or absent** → the module really is variable throughout,
  the `docs/10` §6 logic (conservation picks the question, not the answer) applies in full, and
  E93A-style UC5-specific residues stay the better bet.
- **The invariant positions cluster on one face in the AF3 model** → that face is the candidate
  trigger-binding surface and the XL-MS design in `docs/08` P2 should be aimed at it.
