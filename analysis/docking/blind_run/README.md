# Blind ab-initio run (`1_aRES_T7`)

`job_params_blind_v2.json` is the parameter file to submit. Keep it — HADDOCK deletes
results from the server after two weeks, so this plus the output is the only
reproducibility record.

## The sampling funnel

The validated ab-initio recipe is **10000 / 400 / 400** (it0 / it1 / water), used by the
Bonvin lab's own [CASP-CAPRI T70 ab-initio tutorial](https://www.bonvinlab.org/education/HADDOCK-CASP-CAPRI-T70/)
— blind docking with no prior information, i.e. exactly this case. The default
1000 / 200 / 200 is for *restrained* docking.

| Stage | Default | Form export | **Here** | Note |
|---|---|---|---|---|
| `structures_0` (it0, rigid body) | 1000 | 15000 | **15000** | more surface coverage; fine above the recommended 10000 |
| `structures_1` (it1, semi-flexible) | 200 | 200 | **400** | |
| `waterrefine` (water) | 200 | 50 | **400** | |
| `anastruc_1` (analysed/clustered) | 200 | 50 | **400** | |

**Why the narrow funnel mattered.** At 15000 / 200 / 50, only 1.3% of rigid-body solutions
reached flexible refinement and 50 structures reached clustering. In ab-initio mode the it0
score is weakly discriminating — with random AIRs there is no restraint energy to guide the
ranking — so an aggressive cut at it0 risks discarding the correct solution before flexible
refinement, which is the stage where it would have improved. And 50 structures with
`clust_size: 4` gives too few clusters to judge convergence by population. Raising it0 to
15000 only pays off if the funnel below it is wide enough to use the extra sampling.

**If the queue is too long**, reduce `structures_0` to 10000 rather than narrowing the lower
stages again — 10000 / 400 / 400 is the validated combination. it1 is the expensive stage
(flexible MD per model), so 200→400 roughly doubles the costly part.

## Other settings

Unchanged and fine: `ranair: true` (blind), FCC clustering 0.6 / min 4, `filter_buried` 15%,
aRES trimmed to 1–447 (448–459 is disordered, mean pLDDT 53.9).

`iniseed: 917`. For a replicate run, change only this — two seeds tell you whether a cluster
is reproducible or an artifact of one trajectory.

Submit via **Submit file** in the HADDOCK 2.4 navbar, not the web form, or the form will
overwrite these values.

## Analysis when it finishes

```bash
python3 ../../scripts/interface_contacts.py <run>/structures/it1/water/*.pdb \
    --chains B A --domain A:1-120:NTD1 --domain A:126-447:RES \
    --roi B:100-125:loop --tsv ../../../results/blind_contacts.tsv
```

Read every cluster, not just cluster 1. Pre-registered readout in
`docs/07_docking_run_review.md` §3 — and a non-converging result is the expected outcome at
this size and is **not** evidence against the hypothesis.
