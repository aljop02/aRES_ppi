# Review of the blind docking run, and a Zn ribbon in aRES

Review of the HADDOCK 2.4 parameter file `job_params.json` (run `ares_T7_a`), plus a structural
finding that came out of inspecting the input models.

---

## 1. 🔬 aRES has a Zn(Cys)₄ site at its N-terminus

Found by inspecting the aRES AlphaFold model supplied as partner 1. Four cysteines form a tight
tetrahedral cluster:

| | |
|---|---|
| Cysteines | **C8, C11, C31, C34** |
| Sequence spacing | [3, 20, 3] → **C8**xx**C11** … **C31**xx**C34**, a CxxC…CxxC zinc-ribbon pattern |
| All pairwise SG–SG distances | **3.55 – 3.81 Å**, mean **3.67 Å** |
| Deviation from ideal | **0.09 Å** |
| mean pLDDT of the four cysteines | **90.8** |

**The geometry is quantitatively diagnostic.** An ideal tetrahedral Zn(Cys)₄ site with a Zn–S bond
of 2.3 Å gives an S···S edge of 2.3 × √(8/3) = **3.76 Å**. The observed mean of 3.67 Å is 0.09 Å
off that — and AlphaFold is *confident* about it (pLDDT 90.8 across the four residues), so this is
not a low-confidence region being over-read. These are *not* disulfides: an S–S bond is 2.05 Å,
far shorter.

The other three cysteines in aRES (C78, C96, C249) are isolated — no second 4-Cys cluster exists.

**Negative control:** the same test on the T7 gp5 model finds **10 cysteines and no cluster at all**
— the closest SG–SG pair is 7.85 Å, twice the zinc distance. So the method is not simply calling a
site wherever cysteines occur.

This confirms, for aRES specifically, the observation in
[Anantharaman et al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3710099/) that some
RES-fused PF18870 proteins carry an **additional N-terminal Zn ribbon**.

### Why it matters

**The N-terminal block is bimodular.** It resolves the discrepancy noted in `docs/05` §1.1 between
the structural block (1–120) and the Pfam call (49–136):

```
1 ────────── 48 │ 49 ─────────── 120 │ 121-125 │ 126 ──────────── 447 │ 448-459
   Zn ribbon    │   PF18870 (NTD-1)  │  hinge  │   RES side block     │  tail
   C8 C11 C31 C34                                PF08808: 199-363
└──── structural N-block 1-120 (AF3 PAE) ────┘
```

**InterPro's PF18870 call starts at residue 49 — exactly where the Zn ribbon ends.** The two
annotations are not in conflict; they are describing two different modules that the PAE analysis
saw as one compact block.

### Practical consequences

1. **Never mutate C8 / C11 / C31 / C34.** Those are folding residues. A cysteine mutant there
   would abolish defense for trivial reasons and tell you nothing about recognition — exactly the
   failure mode `docs/02` §9 warns about for P5.
2. **A 1–120 construct contains the Zn site.** Zinc availability in expression and purification
   buffers matters; so does avoiding strong chelators. Worth knowing before the pull-down constructs
   in `notes/vault/` are built.
3. **There is no zinc in the docking run.** The model has the geometry but not the metal. Harmless
   for rigid-body sampling, but during flexible refinement that region is unphysical — four
   close, uncoordinated cysteine sulfurs. **If docking contacts land on residues 1–48, discount
   them** until the run is repeated with Zn modelled.
4. **Open question worth keeping:** does the Zn ribbon itself contribute to phage-protein
   recognition? Zn ribbons are frequently protein- or nucleic-acid-interaction modules. If the
   docking interface lands on 1–48 rather than 49–120, that is a different hypothesis, not a
   failure — but see point 3 before believing it.

Reproduce with:

```bash
python3 analysis/scripts/zn_site_check.py <ares_model.pdb>
```

---

## 2. Review of `job_params.json`

### What the file is

Blind ab-initio docking: **`ranair: true`** (random AIRs), no restraint file (`tblfile: None`),
no centre-of-mass restraints (`cmrest: false`), empty active/passive lists.

| | Partner 1 | Partner 2 |
|---|---|---|
| segid | A | B |
| protein | **aRES** (UC5) | **T7 gp5** |
| residues | 459 (1–459, contiguous) | 704 (1–704, contiguous) |
| atoms | 3,740 | 5,622 |
| mean pLDDT | 78.7 | 93.4 |

### ✅ Correct: these are separate monomers, not cofold chains

Both structures are independently centred on the origin — centroids 3.4 Å apart, minimum atom–atom
distance 0.19 Å, i.e. **superimposed, not docked side by side**. That is the signature of two
independent coordinate frames. Chains extracted from a cofold would share one frame and sit
adjacent, with centroids tens of Å apart.

So this is a genuine *unbound* docking test, which is what you want. The overlap is cosmetic —
`randorien: true` randomises orientation and separates the molecules before rigid-body sampling.

### ⚠️ Change before submitting

**1. `structures_0: 1000` → `10000`.** The most important change. 1000 is the default for
*restrained* docking; for ab-initio with random AIRs the Bonvin lab recommends raising rigid-body
sampling to ~10000, scaled to the surface area being searched. Two complete surfaces of 459 and 704
residues is a large search space. This is also what turns it into an overnight job.
([ab-initio docs](https://www.bonvinlab.org/haddock3-user-manual/abinitio_docking.html))

**2. Trim the aRES tail: use 1–447, not 1–459.** Residues 448–459 have mean pLDDT **53.9**, and
MobiDB/your own notes call 438–459 disordered. A floppy 12-residue tail will find spurious contacts
and pollute the clustering. This is construct #2 from `notes/vault/`.

Low-confidence regions, for reference:

| Region | mean pLDDT |
|---|---|
| Zn ribbon 1–48 | 85.6 |
| NTD-1 49–120 | 81.3 |
| **hinge 121–125** | **54.3** |
| RES side 126–447 | 78.4 |
| **flex loop 392–408** | **47.9** |
| **tail 448–459** | **53.9** |

The hinge scoring 54.3 is independent confirmation that it is a real hinge, not a modelling artifact.
The 392–408 loop is internal and cannot be trimmed — just do not interpret contacts there.

**3. `ranair: true` ignores every other restraint.** Per the docs: *"no other type of restraints will
be considered, even if specified in the configuration file."* So this file **cannot** be converted
into the epitope-restrained run by adding a `.tbl`. That run needs `ranair: false` plus active
residues. Worth knowing before generating variants from this JSON.

**4. Vary `iniseed`** (currently 917) between replicate runs. Blind docking is stochastic; two seeds
tell you whether a cluster is reproducible or an artifact of one trajectory.

Settings that are fine as they are: FCC clustering at 0.6 / min 4, `filter_buried: true` at 15%,
`rotate180_0: true`, `anastruc_1: 200`, `waterrefine: 200`.

---

## 3. Blind vs. restrained — the honest answer

Earlier documents argued that the one-sided restrained run is the better design, and then this
review optimised a blind run without restating that. Both positions are defensible, but they need
stating together, because **they answer different questions and the asymmetry matters.**

| | Blind (`ranair`) | One-sided restrained |
|---|---|---|
| What goes in | nothing | the gp5 epitope, from escape-mutant genetics |
| Chance of converging | **low** at this size, with a flexible epitope | moderate |
| If it works | **unimpeachable** — nothing about the answer went in | good, but the gp5 side was assumed |
| If it fails | **tells you nothing** | still localises the aRES side |

The key asymmetry: **a blind run only produces information when it succeeds.** A smeared,
non-converging result is the expected outcome at this size and is not evidence against the
hypothesis — it is just an unsuccessful search. The restrained run, by contrast, yields a usable
answer about the aRES side either way, which is the actual open question.

So the position has not changed, and it was never "do not run blind": from the start the blind run
was worth keeping **as the unbiased control**. What is not a good trade is spending the only
available compute on blind *instead of* restrained, because the most likely outcome is a night with
nothing to show.

**Recommendation: run both.** The blind run is already configured — bump sampling to 10000, trim the
tail, submit it. But do not let it be the only thing running, and decide now what you will conclude
from each outcome:

- Blind converges on gp5 ~100–125 **and** aRES NTD-1 → the strongest computational result available
  in this project. Nothing went in; it came out anyway.
- Blind converges somewhere else, reproducibly across seeds → a real alternative hypothesis worth
  taking seriously.
- Blind does not converge → **no conclusion.** Not a strike against the hypothesis. Move on.

Write that down before looking at the output.

---

## 4. Still cheaper and more decisive than either run

Unchanged from `docs/06` §4 — neither of these depends on docking converging:

1. **AF3 seed replicates.** Is the S122–D104 contact in `notes/vault/` reproducible across 5 seeds,
   or is it one lucky model? Everything structural rests on that one cofold.
2. **Is aRES D104 conserved?** If the equivalent aspartate is present in the T7-responsive homologs
   (UC5, S_34AVR) and absent in GN04826, which does **not** respond to T7 (`docs/06` §1), that is an
   independent evolutionary argument for the interface — and it needs three sequences, not a night
   of compute.
