# Blind ab-initio run (`ares_T7_a`)

`job_params_blind_fixed.json` is the parameter file as actually submitted. Keep it —
HADDOCK results are deleted from the server after two weeks, and this plus the run
output is the only reproducibility record.

Derived from the file exported by the HADDOCK 2.4 form, with three changes
(see `docs/07_docking_run_review.md` §2):

| Setting | Form export | Here | Why |
|---|---|---|---|
| `structures_0` | 1000 | **10000** | 1000 is the default for *restrained* docking; ab-initio with random AIRs needs ~10× that |
| aRES (partner 1) | 1–459 | **1–447** | residues 448–459 have mean pLDDT 53.9 and are disordered; a floppy tail finds spurious contacts |
| `iniseed` | 917 | **1001** | vary between replicates to test cluster reproducibility |

Unchanged: `ranair: true` (blind mode), FCC clustering 0.6 / min 4, `filter_buried`
at 15%, 200 structures for it1 and water refinement.

Submit via **Submit file** in the HADDOCK 2.4 navbar, not the web form.

## Analysis when it finishes

```bash
python3 ../../scripts/interface_contacts.py <run>/structures/it1/water/*.pdb \
    --chains B A --domain A:1-120:NTD1 --domain A:126-447:RES \
    --roi B:100-125:loop --tsv ../../../results/blind_contacts.tsv
```

Read every cluster, not just cluster 1. Pre-registered readout in
`docs/07_docking_run_review.md` §3 — and note that a non-converging result is the
expected outcome at this size and is **not** evidence against the hypothesis.
