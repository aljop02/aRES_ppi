# The aRES metal site: the AF3-with-Zn result, and how to close it out

Follows `docs/07` §1, which found a four-cysteine cluster with metal-site geometry at the aRES
N-terminus from an **apo** model. This document records what the **holo** AF3 run (Zn supplied as a
ligand) added, states plainly what it did *not* add, and gives the plan for turning a prediction
into an established fact.

---

## 1. Where this stands

| | |
|---|---|
| Apo geometry (`docs/07` §1) | 4 Cys, SG···SG 3.55–3.81 Å, 0.09 Å off ideal — *consistent with* a site |
| **Holo AF3 run (this doc)** | **Zn placed on exactly those 4 Cys with correct bond lengths and angles** |
| Independent of structure prediction | **nothing yet** ← this is the gap |
| Measured | nothing |

---

## 2. The AF3-with-Zn result

Run: aRES (UC5, 459 aa) + 1 Zn²⁺, AlphaFold3 server. Files in `analysis/afold/zn_holo/`.

**The Zn went exactly where `docs/07` predicted it would.**

| Measurement | Value | Ideal |
|---|---|---|
| Coordinating residues | **C8, C11, C31, C34** | — |
| Zn–SG distances | 2.04, 2.26, 2.31, 2.49 Å (mean **2.28**) | 2.30 Å |
| S–Zn–S angles | 97.4–119.1°, mean **108.9°** | 109.47° |
| Zn pLDDT | **95.9** | — |
| pLDDT of the four Cys | 93.9–95.2 | — |
| PAE, Zn ↔ each of the four Cys | **1.3–1.9 Å** | — |
| Contact probability, Zn ↔ C8/C11/C31/C34 | **0.98–0.99** | — |
| ipTM (protein ↔ Zn) | **0.97** | — |
| `has_clash` | **0.0** | — |
| Next nearest Cys (C78) | **14.5 Å** — nothing competing | — |

Reproduce:

```bash
python3 analysis/scripts/zn_coordination.py analysis/afold/zn_holo/fold_ares_zn_model_0.cif
```

### What this does establish

1. **It is a site, not a surface sticker.** Four ligands in the first shell, nothing else within
   7.7 Å. An ion dumped on the surface gives one or two contacts and no angular agreement.
2. **The geometry is chemically correct.** Mean Zn–S 2.28 Å against an ideal 2.30 Å; mean S–Zn–S
   108.9° against 109.47°.
3. **It is realistically distorted, not suspiciously perfect.** The angular spread is 21.7°
   (97.4–119.1°). Real Zn ribbons in the PDB show exactly this kind of distortion; a model that
   returned perfect Td symmetry would be more worrying, not less.
4. **AF3 chose this site out of seven cysteines** and the whole 459-residue surface, and is
   confident about it (PAE 1.3–1.9 Å to the four ligands, against a median PAE of 19.1 Å to the
   protein as a whole).

### ⚠️ What this does NOT establish — read before quoting it

**You told AlphaFold there is a zinc.** The run cannot answer "does aRES bind zinc?", because the
answer was an input. It answers a narrower, conditional question: *if* a metal is bound, this is the
site, and the apo cluster in `docs/07` was not an artifact of an empty pocket.

That is a real result, and it is not the result. **Every piece of evidence for this site so far comes
from one method — AlphaFold — run twice.** Steps A–D below exist to fix precisely that.

> **Do not use the apo-vs-holo pLDDT comparison as evidence.** The Zn-ribbon region rises (85.6 →
> 88.2) but the RES side, 64 Å from the Zn, *falls* by 6.6 and the full-length mean falls by 4.5. A
> zinc ion cannot change confidence in a domain 64 Å away, so the two runs differ in seed, MSA or
> model version, and the comparison is confounded. The local rise may well be real; this pair of
> runs cannot show it. If you want that number, re-run apo and holo with the same seed and version
> and compare then.

---

## 3. Tomorrow: ordered by information per hour

Total active time ≈ 2 h, plus unattended compute.

### Step A — the negative control (~30 min, mostly waiting) ★ do this first

**The single most valuable remaining computational step**, because it is what breaks the circularity
above.

Re-run AF3 on **T7 gp5 + 1 Zn²⁺**. gp5 is the control already used in `docs/07` §1: 10 cysteines, no
cluster, closest SG–SG pair 7.85 Å.

Then run the same analysis on the output:

```bash
python3 analysis/scripts/zn_coordination.py <gp5_zn_model>.cif
```

| Outcome | Reading |
|---|---|
| Zn uncoordinated, or 1–2 contacts, no tetrahedral agreement | **The aRES result means what it looks like.** AF3 does not invent sites on demand. |
| Zn placed in a tidy 4-ligand tetrahedral site in gp5 too | **The aRES result is near-worthless** — AF3 places a supplied Zn *somewhere* plausible regardless. Drop the argument and go to Step D. |

Run a second aRES replicate with a different seed at the same time; costs nothing extra in attention.

### Step B — Foldseek on the N-terminal fragment (~20 min)

Structure-against-structure, which is the right class of tool for a structural claim.

- Server: **search.foldseek.com**. Databases: **PDB100** and **AFDB/Swiss-Prot**. Mode: 3Di+AA.
- Submit **two** queries, because 48 residues is at the short end of where Foldseek's E-values stay
  reliable:
  1. residues **1–48** alone (the putative ribbon)
  2. residues **1–120** (the whole PAE block — gives the aligner more to work with; then check
     whether the hit's aligned region actually covers 1–48)
- Fall back to **DALI** (`ekhidna2.biocenter.helsinki.fi/dali/`) for the short fragment if Foldseek
  returns nothing interpretable — it tolerates short queries better.

**Decision rule:** hits annotated as zinc-ribbon / rubredoxin-like / treble-clef, with the aRES
cysteines aligning onto the template's metal-binding cysteines.

**The step most people skip, and the one that pays:** if a hit is a *holo* structure, superpose it on
the aRES model and check whether the template's metal lands inside the cavity defined by the four SG
atoms. That converts "similar fold" into "the metal goes in the same place" — an independent
structural confirmation of Step 2 that owes nothing to AlphaFold.

### Step C — HHpred (~30 min, mostly waiting)

Profile-profile, far more sensitive than the RPS-BLAST behind CDD. This is how zinc ribbons normally
get found when CDD misses them (see `docs/07` §1 for why CDD missed it).

- Server: **toolkit.tuebingen.mpg.de/tools/hhpred**
- Databases: **PDB_mmCIF70**, **Pfam-A**, **SCOPe70**
- MSA generation: **HHblits vs UniRef30**, 3 iterations
- Submit **both** full-length 1–459 and **1–60 alone** — boundary and composition effects are exactly
  what is in question here.

**Decision rule: probability > 80 %** to a zinc-ribbon or rubredoxin-like family over residues 1–48.

That threshold is not arbitrary — it is the criterion the Osterman paper itself used to call aRES
systems (`docs/01` §"Detection of putative aRES systems": hits to both PF08808 and PF18870 at
probability > 80 %). Using the same pipeline and the same cutoff makes the result directly quotable
in a methods section alongside theirs.

### Step D — cysteine conservation (~45 min) ★ the one that owes nothing to structure prediction

Cheapest decisive test in the list, and the only one here whose answer is independent of AlphaFold,
Foldseek and HHpred alike.

1. Pull aRES sequences for **UC5**, **GN04826**, **S_34AVR** (and anything else from Tables S2/S3
   once acquired — still open in the README status list).
2. Align: `mafft --localpair --maxiterate 1000` (L-INS-i).
3. Read off the four positions equivalent to **C8, C11, C31, C34** (UC5 numbering).

| Outcome | Reading |
|---|---|
| All four invariant across divergent homologs | **Strong.** Four invariant cysteines in a CxxC…CxxC arrangement across homologs with different phage ranges is hard to explain any other way. |
| Present in some, absent in others | The ribbon is a **modular accessory**, not core to the fold — and then *which* homologs have it, mapped against the Fig. 1B phage range, becomes an interesting question in its own right. |
| Absent in GN04826 specifically | Note it loudly: GN04826 is the one that does **not** defend against T7 (`docs/06` §1). That would be a different hypothesis, not a footnote. |

Run this in the same alignment as the **D104 conservation test** already queued in `docs/07` §4.2 —
same sequences, same alignment, two readouts. There is no reason to do them separately.

### Step E — the measurement (order reagents tomorrow; run when protein exists)

None of A–D measures metal. These do, roughly in increasing order of effort:

| Method | What it gives |
|---|---|
| **Thermal shift (DSF) ± EDTA** | Cheapest. A structural metal gives a sharp Tm drop on chelation. Needs only the DSF instrument and purified 1–120. |
| **PAR/PMPS colorimetric assay** | Zn stoichiometry per monomer, on a plate reader. Standard, cheap, quantitative. |
| **Native / intact-mass MS** | A clean +65 Da on the intact mass. Elegant and unambiguous if you have access. |
| **ICP-MS** | The reference method, and it identifies the metal rather than assuming it — which is what settles the Zn-vs-Fe question in `docs/07`. |

Practical note that matters **before** any of this: the 1–120 construct contains the site, so
**keep EDTA out of the purification buffers** and consider 10–50 µM ZnCl₂ in the lysis buffer. A
protein stripped during purification will give a negative answer for a trivial reason.

---

## 4. Register the interpretation before you look

Consistent with the rule in `docs/07` §3 for the docking run. Write these down tomorrow morning:

- **A says "no site in gp5" + B or C returns a zinc-ribbon fold** → the site is established as well
  as computation can establish it. Write it up and move to Step E.
- **A says "no site in gp5" but B and C both come back empty** → the metal site is probably real but
  the *fold* is not a canonical ribbon. Keep the site, **stop calling it a "Zn ribbon"**, and call it
  a tetrathiolate metal site. The practical consequences in `docs/07` §1 are unaffected.
- **A places a tidy site in gp5 too** → the AF3 evidence collapses. Step D alone decides.
- **D says the cysteines are not conserved** → this is a UC5 idiosyncrasy, not a feature of aRES
  systems. Interesting, but it stops being load-bearing for the project.

---

## 5. What actually changes downstream

This is not curiosity-driven; three live decisions depend on the answer.

1. **Construct boundaries.** `docs/05` §1.1 lists N-block **1–120** as construct #3. If the site is
   real, 1–120 is right and a 49–120 construct — which the InterPro boundary would suggest — would
   be missing half its metal ligands and would almost certainly misfold. The CDD boundary of 29
   (`docs/07` §1) would be worse still: it cuts between C11 and C31.
2. **Buffers.** No EDTA; consider supplementing Zn. Applies to the pull-downs in `notes/vault/` as
   well as to anything structural.
3. **Docking interpretation.** `docs/07` §1 says to discount contacts landing on 1–48 because the apo
   model has four uncoordinated sulfurs there. With a validated holo model, **re-run the docking with
   Zn present** and that caveat goes away — and the open question in `docs/07` §1.4 (does the ribbon
   itself contribute to phage-protein recognition?) becomes answerable instead of unanswerable.

And the one that is not downstream of anything: **never mutate C8/C11/C31/C34 as a test of
recognition.** That holds already and does not wait on any of this.
