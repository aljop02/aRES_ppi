# Three aRES orthologs cofolded with Zn

Run record and provenance:
[`analysis/afold/cofold_zn_three_species/SOURCE.md`](../analysis/afold/cofold_zn_three_species/SOURCE.md).
Machine-readable output: [`results/zn_cofold_sites.tsv`](../results/zn_cofold_sites.tsv),
[`results/zn_cofold_identity_profile.tsv`](../results/zn_cofold_identity_profile.tsv).

**Verdict. The Zn(Cys)₄ ribbon that [`docs/07`](07_docking_run_review.md) §1 inferred from UC5
cysteine geometry alone is confirmed with the metal actually present, and it is the same site in
all three orthologs — 15/15 models, no exceptions.** Beyond that, these runs are mostly *negative*
work on the project's own claims: they qualify the charge-split finding as UC5-specific, show D104
is not conserved, and retract an apparent domain rearrangement as a PAE artifact. One result is new
and cheap to act on: NTD-1's acidic character is a family feature while its sequence is the
*least* conserved part of the protein.

These are single-chain predictions with one ion. They say nothing directly about phage-protein
recognition — no trigger was present in any job.

---

## 0. What was run

| job | label | length | Zn | models | ranking | pTM | chain pTM | ipTM | frac. disordered |
|---|---|---|---|---|---|---|---|---|---|
| `fold_ares` | aRES (UC5) | 459 aa | 1 | 5 | 0.960 | 0.60–0.61 | 0.60–0.61 | 0.97 | 0.12–0.13 |
| `fold_ares_echengduensis` | aRES_Echengduensis | 459 aa | 1 | 5 | 0.960 | 0.58 | 0.58 | 0.97 | 0.13 |
| `fold_kleb_ares` | Kleb_aRES | 443 aa | 1 | 5 | 0.960–0.970 | 0.70–0.71 | 0.70–0.71 | 0.98 | 0.08–0.09 |

The `aRES` job is byte-identical to `data/sequences/ares_uc5.fasta`, so **UC5 numbering applies
with no offset** and every residue label below means what it means elsewhere in the repo.

**Do not read the ipTM of 0.97–0.98 as an interface score.** With a single ion as the second
"chain", ipTM says the model places the ion confidently relative to the protein. It does not
validate the protein, whose own pTM is 0.58–0.71 in the same models. Both numbers are honest and
they are about different things; the gap between them is the whole reason this section lists them
separately.

⚠️ The two ortholog accessions were not recorded with the submission — see SOURCE.md. Nothing here
depends on the species names being right, but they must be pinned down before this leaves the repo.

---

## 1. 🔬 The Zn site is real, and it is the same site in all three

| | UC5 | E.chengduensis | Klebsiella |
|---|---|---|---|
| coordinating residues | **C8, C11, C31, C34** | **C8, C11, C31, C34** | **C8, C11, C31, C34** |
| models with that shell | **5/5** | **5/5** | **5/5** |
| mean Zn–SG | **2.27 Å** | 2.28 Å | **2.22 Å** |
| Zn–SG range | 2.04–2.53 Å | 1.71–3.13 Å | 2.06–2.39 Å |
| tetrahedral RMSD (mean) | 7.2° | 9.2° | **3.2°** |
| tetrahedral RMSD (range) | 4.5–11.6° | 4.7–18.4° | 2.2–5.8° |
| Zn pLDDT | 95.6–96.0 | 96.4–96.6 | 96.8–97.0 |
| Zn ↔ ribbon PAE | 3.6 Å | 3.1 Å | 3.2 Å |

Reference Zn–S(Cys) is **2.30 Å**; ideal tetrahedral angles are 109.47°. The observed means sit
0.02–0.08 Å off the bond length, and the Klebsiella site is within 3.2° of ideal geometry.

This closes the loop on `docs/07` §1, which predicted the site from SG···SG distances (mean 3.67 Å
against an ideal 3.76 Å) in a model that contained **no metal**. With Zn supplied, the metal goes
exactly where that geometry said it would, in every model of every ortholog.

**Why the geometry is quoted rather than the ipTM.** AF3 will place a supplied ion *somewhere*
regardless — a confident-looking ion is the default, not evidence. A site is only credible when the
geometry is right, the ligands are confident, the shell recurs across models, and the PAE pins the
ion to the ligands. All four are reported above, so none of them is being assumed:
`analysis/scripts/af3_zn_site.py` prints all four by construction.

**QC, one bad model of fifteen.** E.chengduensis model 2 is a distorted site: Zn–SG(C8) at 1.71 Å
(impossibly short) and C31 pushed out to 3.13 Å, tetrahedral RMSD 18.4°. It is the *same four
cysteines*, just badly refined, and it is the sole reason that ortholog's range looks wide. The
other 14 models are clean. Worth knowing because it is also a reminder that one AF3 model is not a
measurement.

### Consequences

1. **`docs/07`'s "never mutate C8/C11/C31/C34" now holds family-wide**, not just for UC5. These are
   folding residues in every ortholog examined.
2. **Zinc in every buffer, for every ortholog construct.** The EMBL Hamburg constructs in
   [`docs/08`](08_embl_hamburg_structural_plan.md) all contain residues 1–48, so this is a
   purification constraint for the whole family, not a UC5 quirk. Avoid strong chelators.
3. **The ribbon and NTD-1 are one module.** Ribbon ↔ NTD-1 PAE is 5.0–6.4 Å across the orthologs —
   determined, if only in outline. The bimodular-N-block picture in `docs/07` §1 survives contact
   with two more sequences.

---

## 2. The ribbon is the most conserved module; NTD-1 is the least

Percent identity to UC5, over a global pairwise alignment (BLOSUM62, affine gaps), by region:

| ortholog | overall | Zn-ribbon 1–48 | NTD-1 49–120 | RES side 126–447 | tail 448–459 |
|---|---|---|---|---|---|
| E.chengduensis | 71.9% | 70.8% | 69.4% | 73.6% | 50.0% |
| Klebsiella | 44.5% | 47.9% | **33.8%** | 47.4% | 33.3% |

In the distant ortholog, **NTD-1 is 13.6 points less conserved than the RES domain that follows it
and 14.1 points less than the ribbon that precedes it** (NTD-1 vs RES: z = −2.08, p = 0.037,
two-proportion test). Three independent signs point the same way:

- **Indels avoid the ribbon entirely.** Klebsiella's only two N-block indels are a +1 insertion at
  ref~63 and a −1 deletion at ref~108 — **both inside NTD-1**, and they cancel, so the 49–120
  boundaries hold in all three frames. The ribbon takes no indels in either ortholog.
- **Invariant positions:** 20/48 (42%) of the ribbon is identical in all three, against 21/72 (29%)
  of NTD-1.
- The four ribbon cysteines are at **identical positions** (8, 11, 31, 34) in all three despite
  Klebsiella diverging to 44.5% overall.

**This is consistent with NTD-1 being a specificity-determining module and not much else.** A
sensor that must track different phage proteins in different hosts should diverge faster than the
catalytic domain it activates and the structural ribbon it sits on — which is what the table shows.

**How much weight this carries: not much yet, and the statistics are weak on purpose to say so.**
p = 0.037 is one pairwise comparison, uncorrected, over n = 71 positions, with no phylogenetic
correction and no correction for the several comparisons in the table. Two of the three sequences
are 71.9% identical to each other, so this is effectively **one** independent comparison. It is
suggestive and it is cheap to make rigorous — see §7.

---

## 3. NTD-1's acidic character is conserved; the charge *split* is not

Lys+Arg against Asp+Glu, per region, per ortholog (His never counted — its charge depends on
environment, not composition):

| region | UC5 | E.chengduensis | Klebsiella |
|---|---|---|---|
| Zn ribbon 1–48 | **net +11** (+22.9/100) | net +3 (+6.2/100) | net +2 (+4.2/100) |
| NTD-1 49–120 | **net −16** (−22.2/100) | **net −14** (−19.4/100) | **net −18** (−25.0/100) |
| RES side 126–447 | net −3 (−0.9/100) | net +8 (+2.5/100) | net −8 (−2.5/100) |

**NTD-1 is strongly acidic in all three, and the composition is conserved far beyond the sequence:
zero lysines and exactly eleven aspartates in all three**, in a 72-residue window that is only
33.8% identical in Klebsiella. Seven or eight glutamates each. A module can lose two-thirds of its
sequence identity and keep its charge to within four units — that is a property under selection.

**But the opposite-charge split is a UC5 feature.** The ribbon's +11 drops to +3 and +2 in the
orthologs, so the equal-and-opposite pairing of the two N-terminal modules — the finding in the
most recent commit, and the basis of `charge_composition.py`'s docstring — **does not generalise**.
The acidic NTD-1 is the family-level fact; the dipole is not.

Two consequences, pulling in opposite directions:

- **It strengthens `docs/09` §5 where that rests on NTD-1.** The electrostatic-complementarity
  reading of cluster 6 — acidic NTD-1 against the basic gp5 loop — is based on the one property
  these runs show to be conserved. That experiment gets a better prior.
- **It weakens `docs/09`'s nucleic-acid-partner argument for the ribbon.** That argument ("the
  charge model says the ribbon's partner is nucleic acid", `docs/09` next-steps item 4) needs the
  ribbon to be basic. In two of three orthologs it is near-neutral. The DNA-containing re-dock is
  still worth running — gp5 *is* a DNA polymerase, which is reason enough — but the charge
  argument for it should not be quoted as family-level support.

---

## 4. ❌ Withdrawn: the apparent domain rearrangement

Superposing each ortholog on UC5 one domain at a time gives low RMSDs — Klebsiella's ribbon 1.31 Å,
NTD-1 1.91 Å, RES 2.79 Å — while the **whole chain gives 8.05 Å**. Fitting on one domain and
measuring another in that frame makes it look dramatic:

| fit on | measured | Klebsiella | E.chengduensis |
|---|---|---|---|
| Zn-ribbon | RES | 29.2 Å | 6.2 Å |
| NTD-1 | RES | 26.2 Å | 6.1 Å |
| RES | Zn-ribbon | 20.7 Å | 2.5 Å |

Read naively that is a different hinge angle — the modules conserved, their arrangement not. **It
is an artifact, and the block PAE says so plainly:**

| mean PAE (Å), 5 models | ribbon ↔ NTD-1 | N-block ↔ RES core | RES core internal | tail ↔ anything |
|---|---|---|---|---|
| UC5 | 5.8 | **18.7** | 8.9 | 28–31 |
| E.chengduensis | 5.0 | **22.8** | 12.9 | 29–31 |
| Klebsiella | 6.4 | **12.5** | 8.0 | 21–22 |

**AF3 never determined the N-block-to-RES orientation in any of the three jobs.** At 12.5–22.8 Å
inter-block PAE the relative placement is a guess, so a cross-frame RMSD between two guesses
measures the guessing, not a structural difference. The 20–29 Å numbers above are withdrawn, and
the whole-chain RMSDs (2.88 Å and 8.05 Å) are uninterpretable for the same reason.

This is the check that has to run before *any* inter-domain claim from these models, which is why
`analysis/scripts/af3_pae_blocks.py` exists and why it prints a verdict line per block pair. Two
further things it shows:

- **The tail is as bad as `docs/07` §2 said.** Residues 448–459: pLDDT 49.9–51.4, PAE 28–31 Å
  against everything. Trim it; the advice stands.
- **`postRES` 364–447 is also poorly placed** (15–26 Å PAE to the N-block). The well-determined
  part of the RES side is the PF08808 core, 199–363.

---

## 5. D104 is not conserved

| UC5 | E.chengduensis | Klebsiella |
|---|---|---|
| **D104** (pLDDT 63) | **D104** (pLDDT 63) | **G105** (pLDDT 63) |

The D104–S122 hydrogen bond is the structural basis of the two-escape-mutation story in
[`docs/05`](05_status_after_vault_integration.md) §2 and the README's point 1. It is now resting on
three weak legs at once:

1. **one** cofold model, with no replicate seeds;
2. a residue at **pLDDT 63** in every ortholog — the lowest-confidence position of any discussed in
   this document;
3. **no independent support** from the blind run, which never put D104 within 16 Å of S122
   ([`docs/09`](09_blind_zn_run_results.md) §4);
4. and now a position that is **not conserved** — glycine in Klebsiella, which cannot donate or
   accept that bond at all.

**This does not refute the contact.** Klebsiella aRES may have a different trigger entirely —
`docs/06` §2 already establishes that at least one aRES trigger is not a polymerase, so there is no
reason its NTD-1 should contact a T7 gp5 loop. A specificity module *should* vary at its contact
residues. But the contact can no longer be described as a general recognition mechanism, and the
**AF3 seed replicates stay the load-bearing next step** they became in `docs/09`.

---

## 6. E89 and E93 test different hypotheses — run both, for different reasons

`docs/09` §5 proposed aRES **E89A / E93A** plus gp5 R111A. These runs separate the two:

| UC5 | E.chengduensis | Klebsiella | reading |
|---|---|---|---|
| **E89** | E89 (identical) | **D90** (acidic) | acidic in all three |
| **E93** | **S93** | **N94** | charge lost in both orthologs |
| D90 | D90 | D91 | **invariant in all three** |

The tempting conclusion is "prioritise E89A because it is conserved". **That is the wrong filter
here, and worth stating explicitly, because §2 argues NTD-1 is a specificity module.** If NTD-1
really does track different phage proteins in different hosts, then its trigger-contacting residues
are the ones *expected to vary*, and conservation selects against exactly the residues you want.
So:

- **E89A** (acidic in all three) tests a **general structural or activation** role. A null result
  here is informative — it rules out a conserved acidic position.
- **E93A** (charge lost in both orthologs) tests a **T7/UC5-specific recognition** role. Under the
  sensor hypothesis this is the *more* interesting mutant, not the weaker one.
- **D90A** is worth adding as a third: invariant across all three, immediately adjacent to E89, and
  it would distinguish "any acidic residue in this patch" from "this specific one".

All three are single site-directed mutants in the plasmid-based 4G/4H assays, so the marginal cost
of the third is near zero. None of them is licensed by conservation alone — conservation here picks
the *question*, not the answer.

---

## 7. What these runs do and do not license

**Licensed:**
- A Zn(Cys)₄ ribbon at C8/C11/C31/C34 in all three orthologs, with correct metal geometry (§1).
- Ribbon + NTD-1 as one structural module, PAE 5.0–6.4 Å (§1).
- NTD-1 as a strongly acidic, zero-lysine module in all three (§3).
- NTD-1 as the fastest-diverging of the three modules — suggestively, at p = 0.037 on effectively
  one independent comparison (§2).

**Not licensed:**
- Anything about the N-block-to-RES architecture, hinge angle, or whole-chain shape (§4).
- Anything about the tail, or about 364–447 (§4).
- The D104–S122 contact, which these runs weaken rather than support (§5).
- Any claim about phage-protein recognition: **no trigger was in any of these jobs.**

**A free folding-control set for the chimera experiment.** The README's P2 — swapping NTD-1 between
orthologs — now has two donors with known properties. **Klebsiella NTD-1 is the informative donor**
(33.8% identity, two indels, the real test of modularity) **and the one most likely to misfold**;
**E.chengduensis NTD-1 at 69.4% is the folding control** that should swap cleanly. Running both
distinguishes "the swap failed because NTD-1 is not modular" from "the swap failed because the
construct is broken", which a single chimera cannot.

**And a construct note:** Klebsiella aRES is **naturally the trimmed construct**. Its 13-residue
deletion at ref~438 removes precisely the segment that `docs/07` §2 says to cut from UC5, and its
C-terminus (426–443) sits at pLDDT 70 against UC5's tail at 50.6. If a PF18870-plus-RES protein is
wanted that needs no truncation design, that is the one.

---

## 8. Next steps, cheapest first

1. **Record the two accessions** in SOURCE.md and add both sequences to `data/sequences/`. They are
   marked VERIFY until then, and §2's whole comparison is unciteable without them.
2. **E89A, E93A and D90A** through the 4G/4H assays, framed as in §6 — they test different
   hypotheses, so run them together and do not treat E93A as the weak one.
3. **AF3 seed replicates of the aRES×gp5 cofold.** Unchanged from `docs/09`, and §5 raises the
   stakes: D104 is now both low-confidence and non-conserved.
4. **Zn in all buffers** for every ortholog construct (§1).
5. **Run `af3_pae_blocks.py` before any inter-domain claim** from an AF3 model, including the
   existing cofold. §4 is what happens when that check is skipped.
6. **Use both orthologs in the chimera panel** (§7), Klebsiella as the informative donor and
   E.chengduensis as the folding control.
7. **Ten to twenty PF18870+RES orthologs** from the Pfam alignment, same per-region identity
   analysis with phylogenetic correction. This is the one that turns §2 from suggestive into
   evidence, and it is pure compute — no wet work and no new predictions needed.

### Registered in advance, for the ortholog panel in item 7

- **NTD-1 identity significantly below both flanking modules**, with phylogenetic correction → the
  modular-sensor hypothesis gains real comparative support, and NTD-1's variable positions become
  the primary mutagenesis targets.
- **NTD-1 as conserved as the RES domain** → the hypothesis loses its comparative leg. It would not
  be refuted (a sensor can be conserved if the triggers are), but §2 could not be cited for it, and
  the case would rest entirely on the functional assays.
- **The ribbon diverging as fast as NTD-1** → the bimodular reading in `docs/07` §1 is wrong, and
  the N-block should be treated as one evolutionary unit rather than two.
