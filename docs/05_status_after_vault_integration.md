# Status after integrating the existing work (Oct 12, 2026)

Your Obsidian vault is now in `notes/vault/`. This document reconciles it with `docs/01`–`04`,
corrects what those got wrong, and says what I think the decisive next steps are.

**Headline: your cofold already supports the hypothesis, and the two escape mutations fit together
into one coherent mechanism rather than two competing ones.** There is also a measurement
inconsistency across the three PyMOL panels that should be fixed before the numbers are used.

---

## 1. What the vault establishes that `docs/01`–`04` did not have

### 1.1 Real domain boundaries — and the structural module ≠ the Pfam call

From your InterPro/PAE analysis (`notes/vault/Central question...md`,
`notes/vault/figures/uc5_ares_pae_annotated.png`):

| Source | Boundaries |
|---|---|
| InterPro **PF18870** | **49–136** |
| InterPro **PF08808** RES | **199–363** |
| MobiDB disordered | 438–459 |
| **AF3 PAE structural module** | **N-term block 1–120** · hinge 120–125 · RES-side block 126–447 |

> **Correction (and it reverses what this document first said).** The RES Pfam accession is
> **PF08808**, which is what your note and your PAE figure already had. The paper is internally
> inconsistent: Fig. 1B labels the domain `PF08808 RES`, while the STAR Methods text repeatedly
> writes `PF00808` — a histone-like family, not RES. An HMM search built on the Methods number
> returns the wrong family and would silently corrupt the homolog set for `docs/03` T4. See
> `docs/08` §"Two corrections".

**The important point: the structural module (1–120) and the Pfam annotation (49–136) disagree at
both ends.** Residues 1–48 belong structurally to the N-terminal block but fall outside PF18870,
and PF18870's C-terminal 16 residues (121–136) cross the hinge onto the RES side. Anyone who builds
constructs from the Pfam coordinates alone would cut through the middle of a folded module and
through the hinge.

**→ My `docs/03` proposed NTD-1 constructs as "1–150". That was a guess and it is wrong. Use your
PAE-derived boundaries (1–120 or 1–124 / 126–447).** `docs/03` has been corrected.

Your construct list is the right one and I would not change it:

1. Full length 1–459
2. Tail-trimmed 1–447 (conservative; 438–447 still packs on the RES side)
3. N-block 1–120 (or 1–124)
4. RES side 126–447 (retains R182, R184)

### 1.2 Your cofold already places the escape residue against NTD-1

`notes/vault/figures/Pasted image 20261010161850.png`:

```
gp5 Ser122 OG  ——  3.062 Å  ——  aRES Asp104 OD2
    (chain A)                        (chain B)
```

**aRES Asp104 lies inside the N-terminal block (1–120) and inside PF18870 (49–136).**

So your own cofold independently produces exactly what `docs/02` predicted from the escape-mutation
genetics (P1) and from domain logic (H1): the residue whose mutation abolishes activation sits in
direct hydrogen-bonding contact with the uncharacterized N-terminal domain. Two independent lines —
phage genetics and structure prediction — converging on the same residue pair is the strongest
result in the project so far.

Chain assignment is confirmed internally: the WT panel shows `OD2/106 ASP/A — N/121 GLY/A`, i.e.
D106 and G121 both in chain A, which is gp5. Your MD plan's step 5 (confirm chain assignment) is
still worth doing against the actual files, but the labels are self-consistent.

### 1.3 1T7P does not contain the loop

From your MD plan: **1T7P carries a deletion of residues 118–123**, so Gly121 and Ser122 are absent,
and residues 120–150 are a T7-specific insertion relative to Klenow — flexible and poorly templated.

**→ `docs/04` listed 1T7P as "the docking input". That was wrong and is now corrected.** 1T7P
validates the rest of the fold; the loop that the entire hypothesis turns on exists only in AF
models. This is a genuine weakness of the structural evidence and should be stated plainly in any
write-up: *the interface residue is in a region with no experimental structure.*

It also reframes my "insertion1" lead from `docs/02` §3. The literature element I found
(helix-loop-helix, ~102–122, exonuclease domain, proofreading) and your loop are the same region,
but the useful description is yours: a T7-specific, poorly ordered insertion that the crystal
structure could not resolve. The epitope claim survives; the "well-defined helix-loop-helix"
framing does not.

### 1.4 A specificity statement I could not verify

Your central-question note says aRES is activated by T7 in *E. coli* UC5 and *K. pneumoniae*
"but notably not from T7 of *E. chengduensis*".

I could not confirm this from the paper's main text — it would have to come from the Fig. 1B
heatmap. **If GN04826 aRES genuinely does not defend against T7 while UC5 and S_34AVR do, that is
the single most valuable fact in the project**, because it hands you a matched natural pair for the
domain-swap experiment: one aRES that responds to T7 and one that does not, differing in NTD-1.

Please confirm it against Fig. 1B at full resolution before building on it. Flagged in
`docs/04` as an open item.

---

## 2. The two mutations are one mechanism, not two hypotheses

Your MD plan frames H1 (loop conformation) and H2 (direct contact) as alternatives. Looking at the
three panels together, I think they are **sequential steps of a single mechanism**, and that the
data already point that way:

| Residue | Bond seen in your models | Type | Role |
|---|---|---|---|
| **D106** | D106 OD2 → G121 N, 2.713 Å (WT) | **intramolecular**, within gp5 | holds the loop in the conformation that presents S122 |
| **S122** | S122 OG → aRES D104 OD2, 3.062 Å | **intermolecular**, gp5 → aRES | the actual recognition contact |

Read that way:

> **D106 is the scaffold; S122 is the contact.** D106G escapes *indirectly* by collapsing the loop
> conformation that displays S122. S122P escapes *directly* by deleting the hydroxyl donor (and
> adding a proline that rigidifies the backbone).

This explains cleanly why two chemically unrelated mutations, 16 residues apart, produce the same
phenotype — which is otherwise a coincidence that needs explaining. It also makes a sharp
prediction: **D106G should weaken binding, S122P should abolish the specific contact**, so the two
need not behave identically in a pull-down.

Your H3 (neither; the cofold is simply wrong) stays live and must stay live. The loop is in a region
with no experimental structure and AF3 is being asked to model exactly the part that crystallography
could not resolve. Treat the S122–D104 contact as a *lead*, not a result, until a pull-down supports it.

---

## 3. ⚠️ Measurement inconsistency across the three PyMOL panels

The three distance measurements are **not like-for-like**, so the numbers cannot be compared as they
stand:

| Model | Figure | Atom pair measured | Distance |
|---|---|---|---|
| WT | `20261010154010.png` | **OD2**/106 ASP/A → N/121 GLY/A | 2.713 Å |
| S122P | `20261010154432.png` | **CB**/106 ASP/A → N/121 GLY/A | 2.777 Å |
| D106G | `20261010153742.png` | **N**/106 GLY/A → N/121 GLY/A | 6.440 Å |

Three different atoms on residue 106. Specifically:

- **The S122P panel measures CB, not OD2.** That looks like a mis-pick in PyMOL. It also cannot be a
  hydrogen bond — CB is a carbon with no polar hydrogen, and a C···N separation of 2.78 Å would be a
  steric clash, not a contact. Your note next to it ("H-Bond not disrupted?") is reading a number
  that does not describe the H-bond. **Re-measure OD2/106 → N/121 in the S122P model before drawing
  any conclusion.** It is entirely possible the H-bond really is intact in S122P — that would fit the
  two-step model above, where S122P acts on the contact rather than the loop — but the current panel
  does not show it.
- **The D106G panel has no alternative** — glycine has no side chain, so N→N is all there is. But
  that means it cannot be compared to a WT OD2→N distance either.

**Fix: use a side-chain-independent geometric measure across all models**, e.g. **Cα(106)–Cα(121)**
distance and the 100–125 loop backbone RMSD against the WT model. Those are defined for every
variant including D106G, and they measure what you actually care about (loop conformation) rather
than one bond that only exists in some variants. Report the OD2→N H-bond separately, only where Asp
is present (WT, S122P, S122A).

This matters beyond tidiness: it is the kind of thing a reviewer finds, and it is the difference
between "the H-bond is unaffected in S122P" and "we measured the wrong atom".

---

## 4. The experiment that answers your central question cheapest

Your note asks exactly the right question:

> "do these mutants stop binding or only abolish activation?"

That question does not need MD, and it does not need mutant phages. The paper already built the
platform:

- **Fig. 4G — co-expression toxicity.** aRES + T7 DNApol co-expressed on plasmids in *E. coli*
  MG1655 causes growth arrest **with no phage present** (1 mM IPTG + 0.2% arabinose, 25 °C,
  OD₆₀₀ on a plate reader). **Activation readout.**
- **Fig. 4H — pull-down.** His-tagged catalytically dead aRES (R182A or E236A) + T7 DNApol,
  Ni-NTA, SDS-PAGE. **Binding readout.**

Both are plasmid-based, so every gp5 variant is site-directed mutagenesis — **no phage genetics
required.** Run the same variant panel through both assays and the binding/activation question is
answered directly:

| gp5 variant | Purpose | Binding (4H) | Activation (4G) |
|---|---|---|---|
| WT | reference | + | + |
| **S122P** | known escape | ? | − (expected) |
| **S122A** | ⭐ **removes the OG donor without adding proline rigidity** | ? | ? |
| D106G | known escape | ? | − (expected) |
| D106A | charge loss, no added flexibility | ? | ? |
| D106N | isosteric, loses charge, keeps H-bond capacity | ? | ? |
| G121A | tests the acceptor side of the intramolecular bond | ? | ? |

**S122A is the discriminator.** Serine→alanine deletes the hydroxyl that donates to aRES D104 but
adds none of proline's backbone restriction:

- **S122A escapes** → the mechanism is the **direct H-bond** to aRES D104 (your H2). The cofold
  interface is real and you have validated an AF3 prediction by mutagenesis.
- **S122A behaves like WT** → the proline's effect on backbone conformation is what matters (your
  H1), and the specific OG→D104 contact is not essential.

Either outcome is publishable-grade mechanism, and it costs one primer pair.

**The reciprocal test is just as cheap and nobody has done it: aRES D104.** If D104 is the partner,
then **aRES D104A / D104N should lose responsiveness to WT gp5 while retaining NADase activity.**
That is the aRES-side proof of the interface, and it is the first functional test ever aimed at
PF18870. Pair it with R182A/E236A as catalysis-dead controls so you can show loss of *defense* with
retained *catalysis* — see `docs/02` §9 on why that control is non-negotiable.

> If one experiment gets done in the rotation beyond the chimeras, make it **gp5 S122A + aRES D104A
> through the 4G/4H assay pair.**

---

## 5. On the MD plan

It is a good plan and I would not restructure it. Three comments:

1. **Tier 2 (the complex) is more informative than Tier 1 and should not wait on it.** The question
   that matters is whether the S122–D104 contact persists, and whether the WT interface is stable at
   all. If the WT complex interface falls apart in the first 100 ns, you learn that before investing
   in mutant monomer runs.
2. **Add the Cα–Cα and loop-RMSD metrics from §3** as the primary H1 readout, so the analysis does
   not depend on an atom that is absent in some variants.
3. **Your setup step 1 — five AF models with different seeds, keeping models whose 100–125 loop
   differs — is the right instinct and is the single most important step in the plan.** If the loop
   conformation is not reproducible across seeds, the starting structure is a choice rather than a
   prediction, and every downstream number inherits that. Record the spread and report it.

The honest framing for Tier 4 (UC5 vs *E. chengduensis* polymerase): that comparison is only
meaningful once §1.4 is confirmed. Check Fig. 1B first.

---

## 6. Revised priorities

What changed now that the vault is integrated: the cofold exists, the domain boundaries are real,
and the structural lead is already supported. The bottleneck is no longer prediction — it is
**experimental discrimination**.

| Priority | Action | Why |
|---|---|---|
| 1 | **Confirm the Fig. 1B specificity claim** (§1.4) | Gates the chimera design and MD Tier 4; costs one figure lookup |
| 2 | **Re-measure the three panels consistently** (§3) | The current numbers are not comparable |
| 3 | **gp5 variant panel through 4G/4H** (§4) | Answers "binding or activation" directly; no phage genetics |
| 4 | **aRES D104A/D104N** (§4) | First functional test of PF18870; reciprocal proof of the interface |
| 5 | **NTD-1 domain-swap chimeras** (`docs/02` P2) | Still the cleanest test of the central hypothesis; long lead time, start cloning now |
| 6 | **AF3 seed spread on the loop** (MD plan step 1) | Tells you how much to trust the cofold at all |
| 7 | MD Tier 2, then Tier 1/3 | Mechanism, once the above says the interface is worth simulating |
| 8 | Conservation/SDP analysis (`docs/03` T4) | Is aRES D104 conserved? Cheap, independent, and strengthens everything above |

Item 8 deserves a note: **check whether D104 is conserved across aRES homologs.** If the aspartate
at the equivalent position is conserved in T7-responsive aRES proteins and absent or substituted in
non-responders, that is an independent evolutionary argument for the interface, obtained without a
single experiment. It needs the homolog list from Tables S2/S3 — still the top acquisition item in
`docs/04`.

---

## 7. Still outstanding

- Supplementary Tables S2/S3/S4, and Fig. 1B at full resolution (`docs/04`)
- The actual cofold model files — the vault has figures but no `.cif`/`.pdb`. Upload them and
  `analysis/scripts/interface_contacts.py` can give the full interface contact list, the NTD-1 vs
  RES breakdown, and the seed-to-seed reproducibility of the S122–D104 contact in one command.
- The blind HADDOCK run output, for the same reason.
