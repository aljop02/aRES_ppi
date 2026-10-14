# AlphaFold 3 cofolds: three aRES orthologs, each with one Zn

Three AlphaFold Server jobs, downloaded 2026-10-14, analysed in
[`docs/10`](../../../docs/10_ortholog_zn_cofolds.md).

The raw archives (`fold_*.zip`, ~2.8 MB each) sit in this directory and are
**gitignored** (`analysis/afold/**/*.zip`). They are regenerable from the job
requests below, which are reproduced here in full so the runs can be repeated
without the archives.

| archive | job name | chain length | model seed | models |
|---|---|---|---|---|
| `fold_ares.zip` | `aRES` | 459 aa | 1126069797 | 5 |
| `fold_ares_echengduensis.zip` | `aRES_Echengduensis` | 459 aa | 333072434 | 5 |
| `fold_kleb_ares.zip` | `Kleb_aRES` | 443 aa | 54977774 | 5 |

Each job: one protein chain, `useStructureTemplate: true`, plus one `ZN` ion;
`dialect: alphafoldserver`, `version: 3`. Each archive carries paired and unpaired
MSAs, four template hits, and per-model `model_*.cif`, `summary_confidences_*.json`
and `full_data_*.json` (the last holds the PAE matrix).

## Sequence identity of the `aRES` job

The 459-aa `aRES` sequence is **byte-identical** to
[`data/sequences/ares_uc5.fasta`](../../../data/sequences/ares_uc5.fasta)
(*E. coli* UC5). Verified by string comparison, not by eye. Every residue number in
`docs/10` is therefore in the UC5 frame that the rest of the repo uses, with no
offset — and C8/C11/C31/C34, E89, E93 and D104 all refer to the same residues as in
`docs/05`–`docs/09`.

## ⚠️ VERIFY: the two ortholog accessions are unrecorded

The job names are the submitter's labels. `aRES_Echengduensis` and `Kleb_aRES`
presumably mean *Escherichia chengduensis* and a *Klebsiella* sp., but **no
accession, assembly or strain was recorded with the submission**, and neither
sequence is in `data/sequences/`. Nothing in `docs/10` depends on the species
assignment being right — the comparative results hold for "two aRES orthologs, one
close and one distant" however they are named — but the accessions are needed
before any of it is written up outside this repo, and before either sequence is
used to pick a phage panel or a construct.

Per the repo rule in the README, these stay marked VERIFY until the accessions are
recorded here.

## Reproduce the analysis

Extract each archive into its own directory, then:

```bash
python3 analysis/scripts/af3_zn_site.py <dir> --cutoff 3.5 \
    --regions 1-48:Zn-ribbon 49-120:NTD1 126-447:RES 448-459:tail \
    --tsv results/zn_cofold_sites.tsv

python3 analysis/scripts/af3_pae_blocks.py <dir> \
    --blocks 1-48:Zn-ribbon 49-120:NTD1 126-447:RES 448-459:tail

python3 analysis/scripts/homolog_compare.py <uc5.cif> <ech.cif> <kleb.cif> \
    --names UC5 E.chengduensis Klebsiella \
    --regions 1-48:Zn-ribbon 49-120:NTD1 126-447:RES 448-459:tail \
    --map 8 11 31 34 89 93 104 --conserved 1-48:Zn-ribbon 49-120:NTD1 \
    --super 1-48:Zn-ribbon 49-120:NTD1 126-447:RES \
    --profile results/zn_cofold_identity_profile.tsv
```
