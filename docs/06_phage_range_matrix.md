# The phage-range matrix (Fig. 1B) — decoded, and what it breaks

Source: Osterman et al. 2026 Fig. 1B. Values read off the colour scale, so treat adjacent bins
(e.g. 1,000 vs 10,000) as approximate. The qualitative calls — **defends** vs **does not defend** —
are unambiguous and are what the arguments below rest on.

## The matrix

Fold defense: **·** = 1 (no defense) · **+** = 10 · **++** = 1,000 · **+++** = 10,000 · **++++** = >100,000

| aRES from | T2 | T4 | T5 | **T7** | SECphi4 | SECphi18 | SECphi17 | SECphi27 | **λ-vir** |
|---|---|---|---|---|---|---|---|---|---|
| ***E. coli* UC5** | + | +++ | ++++ | **++++** | · | · | · | +++ | **++** |
| ***E. chengduensis* GN04826** | · | ++ | +++ | **·** | · | + | · | ++/+++ | **++** |
| ***K. pneumoniae* S_34AVR** | · | · | · | **++++** | ++ | ++++ | · | · | **·** |

Encodes own DNA polymerase: T2 ✓ (gp43) · T4 ✓ (gp43) · T5 ✓ (own pol, 3′→5′ proofreading exo) ·
T7 ✓ (gp5) · **λ ✗** · SECphi4/17/18/27 unknown.

---

## 1. ✅ Confirmed: GN04826 does not defend against T7

The note in `notes/vault/` was right. **UC5 and S_34AVR both give >100,000-fold against T7;
GN04826 gives 1.** This resolves the open item flagged in `05_status_after_vault_integration.md` §1.4.

### And it is better than just confirmed — UC5 and GN04826 differ at essentially one phage

| | T2 | T4 | T5 | **T7** | SECphi4 | SECphi18 | SECphi17 | SECphi27 | λ-vir |
|---|---|---|---|---|---|---|---|---|---|
| UC5 | + | +++ | ++++ | **++++** | · | · | · | +++ | ++ |
| GN04826 | · | ++ | +++ | **·** | · | + | · | ++/+++ | ++ |
| difference | ~1 log | ~1 log | ~1 log | **~5 logs** | – | ~1 log | – | – | – |

Two homologs with near-superimposable phage ranges, separated by **five orders of magnitude at
exactly one phage — the only phage whose trigger protein is known.**

**This is as close to a controlled natural experiment as you will get, and it makes the domain-swap
experiment (P2) much stronger than I had assumed when I wrote `docs/02`:**

- Build NTD1(UC5)–RES(GN04826) and NTD1(GN04826)–RES(UC5) using the PAE boundaries (1–120 / 126–447).
- Primary readout: **T7 plaque assay**, 5 logs of dynamic range, binary interpretation.
- The other eight phages are **built-in internal controls** in the same assay. T4, T5, SECphi27 and
  λ are defended by *both* parents, so they report on whether the chimeras fold and function at all —
  a chimera that loses everything is a folding failure, not a specificity result. That control comes
  free, and it is exactly what `docs/02` §5 warned would otherwise make a null result uninterpretable.

→ **Use UC5 × GN04826 for the chimeras, not UC5 × S_34AVR.** S_34AVR shares only T7 with UC5, so a
chimera between those two has no internal folding control.

### A supporting detail from the paper

The escape mutants gp5 S122P and D106G escape **both UC5 and S_34AVR** aRES (Fig. 4F). Those two
proteins share **only T7** in their entire range. So two otherwise dissimilar aRES proteins
converge on reading the same loop of T7 gp5 — which is what you would expect if that loop is a
genuine recognition epitope rather than an incidental contact.

---

## 2. ⚠️ λ-vir has no DNA polymerase — so the polymerase model is incomplete

**The observation is correct and it is the sharpest thing in this figure.**

**λ encodes no DNA polymerase.** It replicates with the *host* machinery: the only phage proteins
directly involved are **O** (binds oriλ) and **P** (recruits and remodels the host **DnaB**
helicase). Elongation is done by *E. coli* **DNA polymerase III holoenzyme**, with DnaG primase,
SSB, gyrase and the DnaK/DnaJ chaperones.
([ref](https://pmc.ncbi.nlm.nih.gov/articles/PMC400992), [ref](https://pmc.ncbi.nlm.nih.gov/articles/PMC397919))

Yet **UC5 and GN04826 each give ~1,000-fold defense against λ-vir.** Therefore:

> **At least one aRES trigger is not a phage DNA polymerase.**

This is a materially stronger objection than the authors' own caveat. They wrote that some aRES
defend against T2/T4/T5, "whose DNA polymerases are dissimilar to the T7 polymerase" — but those
phages at least *have* polymerases, so the statement is compatible with NTD-1 reading a divergent
polymerase. **λ has none at all.** No amount of divergence explains it.

So "aRES senses the phage DNA polymerase" is true for T7 and cannot be the general rule.

### What λ could be presenting **[SPECULATION]**

| Candidate | Why plausible |
|---|---|
| **λ P** | Binds and remodels host **DnaB**. The best candidate — see the unifying idea below. |
| λ O | Origin-binding initiator; phage-specific, present early |
| λ Gam | Inhibits host RecBCD — a host-subversion protein, the kind of thing defense systems watch |
| λ Exo / Bet | Red recombination; abundant early proteins |

### A unifying hypothesis worth taking seriously **[SPECULATION]**

Look at what T7 gp5 and λ P have in common:

- **T7 gp5 hijacks a host protein** — *E. coli* **thioredoxin**, as its processivity factor.
- **λ P hijacks a host protein** — *E. coli* **DnaB**, the replicative helicase.

> **NTD-1 may sense phage proteins that capture host replication factors**, rather than polymerases
> as such. That would make it a sensor for a *behaviour* (subversion of host replication machinery)
> rather than for one protein family — which is a far better evolutionary strategy for a defense
> system, because the phage cannot escape by diverging its polymerase sequence.

This connects directly to the thioredoxin point in `02_hypothesis_NTD1_sensor.md` §4: the species
proven to activate aRES in vitro was the **gp5·TrxA complex**, not gp5 alone. If the hypothesis
here is right, that is not an experimental detail — it is the point.

**It is speculation and may well be wrong.** But it is cheap to test, and the test is the assay you
already have.

### ⭐ The experiment

**Co-express individual λ proteins with UC5 aRES and look for growth arrest** — the Fig. 4G format,
exactly as the authors did with T7 DNApol (plasmid-borne, 1 mM IPTG + 0.2% arabinose, 25 °C,
OD₆₀₀ on a plate reader). No phage genetics, no purification.

Order: **P**, then O, then Gam, then Exo/Bet. Use S_34AVR aRES as the negative control — it does
**not** defend against λ (fold defense 1), so a λ protein that activates UC5 but not S_34AVR is
behaving like a genuine specificity-determined trigger rather than a non-specific toxicity artifact.
That control is important: growth arrest on co-expression is easy to get for trivial reasons.

If λ P activates UC5 aRES:
- you have found **the second aRES trigger**, and the first one that is not a polymerase;
- NTD-1 becomes a demonstrated **multi-ligand sensor**, which is the strong form of H2;
- the "host-factor hijacking" model gets its first real support.

Follow up with the pull-down (Fig. 4H) to show direct binding, and co-fold λ P with NTD-1 in AF3.

---

## 3. Two further handles in the matrix

### 3.1 T2 vs T4 — a natural specificity pair that already exists

UC5 gives **10-fold against T2** but **10,000-fold against T4**. T2 and T4 are both T-even phages
with closely related gp43 polymerases. A ~3-log difference between close homologs means the
recognition surface discriminates at fine sequence resolution.

**That is a free experiment.** Align T2 gp43 against T4 gp43, find the differences, and map them
onto the predicted interface and onto the region equivalent to T7's 100–125 loop. If the
discriminating residues cluster there, it is independent support for the epitope — obtained from
sequences that already exist, with no cloning at all.

GN04826 sharpens it further: T2 **1**, T4 **1,000**. Same comparison, cleaner contrast.

### 3.2 SECphi17 is defended by nothing

None of the three aRES touch SECphi17. Useful as a **negative control phage** in any assay, and —
if its proteome is available — as the negative set in a comparative screen: proteins present in
defended phages and absent or divergent in SECphi17.

### 3.3 The three ranges dissect cleanly

| Intersection | Phages | Reading |
|---|---|---|
| UC5 ∩ S_34AVR | **T7 only** | the shared, mapped trigger (gp5) |
| UC5 ∩ GN04826 | T4, T5, SECphi27, λ | shared non-T7 triggers — at least one is not a polymerase |
| GN04826 ∩ S_34AVR | ~nothing (SECphi18 weak vs strong) | the two most dissimilar specificities |

The matrix is doing real work here: it says the T7-gp5 recognition capability is **separable** from
the rest of the range, which is exactly the modularity H2 predicts.

---

## 4. Consequences for the plan

| Change | Where |
|---|---|
| Chimera pair fixed: **UC5 × GN04826**, T7 readout, other 8 phages as folding controls | `docs/02` P2, `docs/03` T6 |
| "aRES senses phage DNA polymerase" demoted from general rule to the T7-specific case | `docs/01` §7, `docs/02` §8.4 |
| **λ protein co-expression screen** added as a top-priority cheap experiment | new, §2 above |
| T2/T4 gp43 comparison added as a zero-cost sequence analysis | new, §3.1 |
| `docs/03` T5 ("the T2/T4/T5 puzzle") is now the **λ puzzle** and is more urgent | `docs/03` T5 |

### Revised top of the priority list

1. **λ P (then O, Gam) co-expression with UC5 aRES** — Fig. 4G format, S_34AVR as negative control
2. **UC5 × GN04826 NTD-1 chimeras**, T7 plaque assay — start cloning now, long lead time
3. **gp5 S122A + aRES D104A** through 4G/4H (`docs/05` §4)
4. T2 vs T4 gp43 alignment — free, do it this week
5. Re-measure the three PyMOL panels consistently (`docs/05` §3)

---

## Sources

- [Osterman et al. 2026, *Cell Host & Microbe* 34:991–999](https://doi.org/10.1016/j.chom.2026.05.004) — Fig. 1B
- [λ replication: O and P proteins, host Pol III holoenzyme](https://pmc.ncbi.nlm.nih.gov/articles/PMC400992)
- [λ P recruits host proteins at oriλ](https://pmc.ncbi.nlm.nih.gov/articles/PMC397919)
- [λO–λP–DnaB complex, replication requirements](https://pmc.ncbi.nlm.nih.gov/articles/PMC3640518)
- [T5 encodes its own DNA polymerase with proofreading exonuclease](https://pmc.ncbi.nlm.nih.gov/articles/PMC256554)
