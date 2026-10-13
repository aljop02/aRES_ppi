# Step-by-step protocol: the one-sided restrained run, with Zn²⁺

Everything needed to submit `2_ares_T7_v2` and know it worked. Written after three runs that each
failed in a different way — two blind runs that did not converge (`docs/09` §2–§5) and one restrained
run that never docked anything (`docs/09` §10).

**The submission package is already built and committed** in
[`analysis/docking/onesided_run_v2/`](../analysis/docking/onesided_run_v2/). You can skip to §3 and
start uploading. §2 is how it was made and how to rebuild it.

---

## 0. What is different this time

| | runs 779823 / 780456 | run 780457 | **this run** |
|---|---|---|---|
| Restraints | random (`ranair`) | **none generated** | **9 explicit AIRs, uploaded as a file** |
| Zn²⁺ in the Zn(Cys)₄ ribbon | absent | absent | **present, restrained** |
| Coordinating cysteines | CYS (neutral) | CYS | **CYF (thiolate)** |
| segid in input PDBs | blank, inferred | blank, inferred | **stamped explicitly** |
| Random restraint removal | n/a | on (50%) | **off** |
| Proof it docked | — | — | **100-structure test job first** |

---

## 1. ⚠️ What the zinc does and does not fix

Worth being straight about this before spending a night of queue time on it.

**What it fixes.** Four uncoordinated sulfurs sitting 3.67 Å apart is not a physical arrangement.
Rigid-body docking does not care, but semi-flexible and water refinement do, and they will distort
that region. With Zn²⁺ present and coordination restrained, the ribbon keeps the geometry AlphaFold
predicted at pLDDT 90.8. It also retires the standing caveat from [`docs/07` §1](07_docking_run_review.md):
*"if docking contacts land on residues 1–48, discount them until the run is repeated with Zn
modelled."* After this run, a 1–48 interface becomes a result you are allowed to believe — which
matters, because "does the Zn ribbon itself contribute to recognition" is a live alternative
hypothesis, not a failure mode.

**What it does not fix — and this corrects what `docs/09` §4 originally said.** The metal-free site
in the previous runs was **not** a "+2 charge hole". HADDOCK's default cysteine is neutral
(protonated SH), so that site was net **neutral**. Modelling it properly gives Zn²⁺ (+2) plus four
CYF thiolates (−1 each) = net **−2**. So adding the zinc makes the ribbon *more* electronegative,
and if anything slightly **more** attractive to gp5's basic DNA-binding cleft, not less.

**The thing that actually fixes the artifact is the AIRs**, which now exist and will drive the
interface. The zinc is for physical realism and for making a 1–48 result interpretable. Do not expect
it to change the scoring behaviour on its own.

---

## 2. Build the inputs (already done — this is the record)

```bash
# aRES: stamp segid A, insert Zn2+ into the C8/C11/C31/C34 site, write its restraints
python3 analysis/scripts/prepare_haddock_input.py \
    <aRES_model>.pdb -o analysis/docking/onesided_run_v2/aRES_Zn_segA.pdb \
    --segid A --chain A --zinc-cys 8,11,31,34 --zinc-resid 500 \
    --tbl analysis/docking/onesided_run_v2/unambig_zn.tbl

# gp5: segid B only
python3 analysis/scripts/prepare_haddock_input.py \
    <gp5_model>.pdb -o analysis/docking/onesided_run_v2/gp5_segB.pdb --segid B --chain B

# the 9 one-sided AIRs
python3 analysis/scripts/make_ambig_tbl.py \
    --from-params analysis/docking/onesided_run/job_params_restrained.json \
    --source-segid B -o analysis/docking/onesided_run_v2/ambig.tbl
```

The two input models were recovered from `partners.N.raw_pdb` inside
`analysis/docking/onesided_run/job_params_restrained.json`, so they are byte-identical to what run
780457 was given.

### What the package contains

| file | contents |
|---|---|
| `aRES_Zn_segA.pdb` | aRES 1–447, segid A, **+ 1 HETATM Zn**, C8/C11/C31/C34 renamed **CYF** |
| `gp5_segB.pdb` | gp5 1–704, segid B |
| `ambig.tbl` | **9** AIRs: each gp5 active residue → all 304 aRES surface residues |
| `unambig_zn.tbl` | **10** restraints: 4 × Zn–SG (2.30 ± 0.10 Å) + 6 × SG–SG tetrahedron edges |

Zinc placement came out at **2.242–2.252 Å** from all four sulfurs — essentially perfectly
equidistant. (Slightly under the 2.30 Å ideal because the observed tetrahedron is a little
compressed: edges 3.55–3.81 Å, mean 3.67 Å, versus 3.76 Å ideal.)

The Zn record, which is the part most likely to be silently mangled:

```
HETATM 3637 ZN+2 ZN2 A 500      40.712  -7.526  17.869  1.00  0.00      A   ZN
```

- atom name `ZN+2` starts in **column 13** (4-character names do not get the usual leading space)
- residue name `ZN2` in columns 18–20
- **chain A — the same chain as the protein.** An ion on its own chain will not move with the
  protein during rigid-body docking and ends up stranded.
- segid `A` in columns 73–76, element `ZN` in 77–78

---

## 3. ⚠️ Verify the ion naming before anything else

This is the one step I could not check for you: this environment's network policy blocks the HADDOCK
and BioExcel sites, so the `ZN2` / `ZN+2` / `CYF` conventions below come from HADDOCK user-forum
threads (consistent across several, but 2.2-era) rather than from the current spec.

**Open the server's supported-ions page and confirm `ZN2` / `ZN+2` is still the expected form.**
Then upload `aRES_Zn_segA.pdb` on its own and check the server accepts it. Three outcomes:

| what happens | what it means |
|---|---|
| accepted, Zn retained | good, continue |
| rejected, "elemental ion without specified charge" | naming is wrong — fix and retry |
| accepted but Zn silently dropped | **the dangerous one.** Check the processed PDB the server returns and count the HETATM records. |

If `CYF` is rejected, rerun the prepare step with `--no-cyf`. You lose the thiolate protonation but
keep the metal and the geometry, which is the larger part of the benefit.

---

## 4. Submit the test job — 100 structures, ten minutes

**Do not submit the full run first.** A test job costs ten minutes and proves the restraints are
live. That is the entire lesson of run 780457, which burned a full 400-structure run proving nothing.

Use the **web submission form (Expert / Guru level)**, not "Submit file → Full parameter file". The
parameter-file route is implicated in the 780457 failure.

**Molecules**

| | molecule 1 | molecule 2 |
|---|---|---|
| file | `aRES_Zn_segA.pdb` | `gp5_segB.pdb` |
| segid | `A` | `B` |
| active residues | **leave empty** | **leave empty** |
| passive residues | **leave empty** | **leave empty** |

Active and passive stay empty *because the restraints are uploaded as files*. Filling both in is how
you end up with two competing definitions.

**Restraints**

- Ambiguous restraints (AIR) file → `ambig.tbl`
- Unambiguous restraints file → `unambig_zn.tbl`
- **"Randomly exclude a fraction of the ambiguous restraints" → NO** (`noecv: false`). The default
  discards 50%; with only 9 AIRs that is reckless.
- Random AIRs (`ranair`) → **off**
- Centre-of-mass restraints (`cmrest`) → **off**

**Sampling (test job only)**

| stage | value |
|---|---|
| `structures_0` (it0, rigid body) | **100** |
| `structures_1` (it1, semi-flexible) | **20** |
| `waterrefine` | **20** |
| `anastruc_1` | **20** |

---

## 5. Check the test job — the four things that matter

Run `2_ares_T7_v2` is only worth scaling up if all four pass.

### 5.1 Did they actually dock? — the primary check

```bash
grep -E 'buried surface area|^REMARK energies' <any cluster pdb>
```

| field | fail | pass |
|---|---|---|
| buried surface area | `-999999` | a positive number, ideally >1000 Å² |
| `REMARK energies` vdw (6th field) | `0` | clearly negative |
| `REMARK energies` elec (7th field) | `0` | clearly negative |

**Buried surface area is the real test**, not Eair. Refining what `docs/09` §10 said: `Eair = 0` on
*separated* partners proves no restraints existed, but on a run that worked the AIRs will be
satisfied and Eair will legitimately be small or zero. So BSA and vdW are what tell you docking
occurred; Eair tells you whether restraints were in force. Check BSA first.

### 5.2 Were the AIRs read?

The run page should report the number of restraints read. Expect **9 ambiguous**, **10 unambiguous**.
If it reports 0 ambiguous, stop — you are about to repeat 780457.

### 5.3 Is the zinc still in its site?

```bash
python3 analysis/scripts/zn_site_check.py <cluster pdb> --regions 1-48:Zn-ribbon
```

Then confirm the Zn–SG distances directly:

```bash
grep -E 'ZN\+2|SG  CYF' <cluster pdb>
```

All four Zn–SG should still be ~2.2–2.4 Å. If the zinc has drifted, the unambiguous restraints are
not being applied — check that `"ZN+2"` kept its double quotes through upload.

### 5.4 Does the gp5 epitope end up at the interface?

```bash
python3 analysis/scripts/interface_contacts.py <run_dir>/*.pdb --chains A B \
    --roi B:102-122:gp5_insertion1
```

The AIRs restrain exactly those residues, so `gp5_insertion1` should be a **large** fraction of the
gp5 interface weight. If it is near the 5–9% the blind runs gave, the restraints are not biting.

---

## 6. The full run

Once §5 passes, resubmit identically with the real sampling:

| stage | value | note |
|---|---|---|
| `structures_0` | **5000** | lower than the blind run — the gp5 side is pinned, so the search space is much smaller |
| `structures_1` | **400** | |
| `waterrefine` | **400** | |
| `anastruc_1` | **400** | |
| `clust_meth` / `clust_cutoff` / `clust_size` | FCC / 0.6 / 4 | |
| `iniseed` | 2002 | |

⚠️ Also check `clust_size` actually takes effect — the blind runs returned clusters of 2–3 models
despite `clust_size: 4` (`docs/09` §3).

---

## 7. Read it out — pre-registered

The readout for this run was registered in
[`onesided_run/README.md`](../analysis/docking/onesided_run/README.md) before any of it ran. Score it
as written, and do not renegotiate it after seeing the output:

| | supported | falsified |
|---|---|---|
| NTD-1 is the binding domain | **≥60%** of aRES interface weight in **1–120** | majority in **126–447** |
| D104 is in the interface | aRES **D104** in the top contacts across clusters | absent from all clusters |

```bash
python3 analysis/scripts/interface_contacts.py <run_dir>/*.pdb --chains A B \
    --domain A:1-48:Zn_ribbon --domain A:49-120:NTD1 \
    --domain A:121-125:hinge --domain A:126-447:RES --domain A:448-459:tail \
    --roi B:102-122:gp5_insertion1 --roi B:258-333:gp5_TBD --roi B:410-704:gp5_polC \
    --tsv results/onesided_v2_contacts.tsv
```

Unlike the blind runs, this one localises the aRES side **whatever the outcome** — that asymmetry is
the whole reason the design exists.

Two things that would still undercut a positive result, and are worth checking before celebrating:

- **The species problem.** gp5 is docked alone, but the activating species in the paper was the
  gp5·TrxA holoenzyme, and TrxA covers the TBD (258–333). An interface touching the TBD is suspect
  until the gp5·TrxA run exists (`docs/02` §4, T1.3 in `docs/03`).
- **The electrostatic artifact** (`docs/09` §4) is reduced by restraints but not eliminated. If the
  aRES side lands on an acidic patch with Eelec in the hundreds, treat it the same way.

---

## 8. If it fails again

| symptom | cause | fix |
|---|---|---|
| BSA `-999999`, vdw/elec exactly 0 | no restraints in force | check the AIR count on the run page; check segids survived upload |
| Zn drifted out of the site | unambiguous restraints not applied | check `"ZN+2"` quoting; check the Zn kept chain A |
| PDB rejected on upload | ion or CYF naming | §3; try `--no-cyf` |
| Docks, but gp5 epitope not at the interface | AIRs read but weak | confirm `noecv: false`; confirm 9 AIRs, not fewer |
| Converges on the Zn ribbon (1–48) | may now be real | the `docs/07` discount no longer applies — but check Eelec before believing it |

**The general rule this project has now paid for three times: a HADDOCK results page full of
clusters, Z-scores and rankings is not evidence that docking happened.** Run 780457 produced 400
ranked clusters in which the two molecules were 270 Å apart. Check BSA before you read anything else.
