# AF3 run: aRES + Zn²⁺ (holo)

AlphaFold3 server run of aRES (UC5, 459 aa) with a single **Zn²⁺** supplied as a ligand. Supplied
Oct 12, 2026. The companion apo run is the model used as HADDOCK partner 1 in
`analysis/docking/blind_run/`.

Purpose: test whether the four-cysteine cluster found in the **apo** model (`docs/07` §1) accepts a
metal in the geometry predicted, and where AF3 puts a Zn when it is given one.

## Files

| File | Tracked | Notes |
|---|---|---|
| `fold_ares_zn_model_0.cif` | no (gitignored) | top-ranked model; Zn is chain B |
| `fold_ares_zn_summary_confidences_0.json` | yes | ipTM 0.97, ptm 0.60, ranking 0.96, `has_clash` 0.0 |
| `fold_ares_zn_full_data_0.json` | no (large, regenerable) | per-token PAE and contact probabilities |

## Result

Zn coordinated by **C8, C11, C31, C34** — Zn–SG 2.04–2.49 Å (mean 2.28 vs ideal 2.30), S–Zn–S mean
108.9° (ideal 109.47°), Zn pLDDT 95.9, PAE to the four ligands 1.3–1.9 Å, no clashes. Next nearest
cysteine is C78 at 14.5 Å.

```bash
python3 ../../scripts/zn_coordination.py fold_ares_zn_model_0.cif
```

## Read this before citing it

**The Zn was an input, not a finding.** This run shows *where* a metal goes if one is present; it
cannot show that aRES binds a metal in vivo. The negative control that makes it meaningful — T7 gp5 +
Zn²⁺ — has **not** been run yet.

Full interpretation, caveats, and the validation plan: **`docs/08_metal_site_validation.md`**.

Do not use this run's pLDDT against the apo run's: the two differ globally (the RES side, 64 Å from
the Zn, drops 6.6) and are not a controlled pair. See `docs/08` §2.
