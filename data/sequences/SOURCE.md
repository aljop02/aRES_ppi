# `data/sequences/` — provenance

Repo rule: no sequence enters this directory without its accession, the date, and the
query that produced it. Anything not pinned down is marked **VERIFY**.

---

## `ares_uc5.fasta` — the reference frame

*Escherichia coli* UC5 aRES, 459 aa. Every residue number in this repository is in this
frame. From **Osterman et al. 2026, *Cell Host & Microbe* 34:991–999**
([10.1016/j.chom.2026.05.004](https://doi.org/10.1016/j.chom.2026.05.004)); see
[`docs/04_data_inventory.md`](../../docs/04_data_inventory.md) for how it was obtained.

**Identified 2026-10-14 — it now has stable identifiers:**

| | |
|---|---|
| UniParc | **UPI00038F43E8** (459 aa) |
| RefSeq | **WP_021575732** — "RES family NAD+ phosphorylase", Enterobacteriaceae, *active* |
| UniProtKB | **none active.** Six TrEMBL entries held this sequence and all are now deleted: `T9RI87`, `V2SFW4`, `A0A136SSX2`, `A0A4V7USR3`, `A0A6D0Z5L7`, `A0A7U8MA60` |
| also in | ~150 *E. coli* genome assemblies (EMBL/EMBL-WGS), incl. `QBZ32856`, `QOX29831`, `WFY50681` |

Method: CRC64 checksum of the sequence → UniParc, then the UniParc sequence pulled back
and compared position by position. **Verified identical, not merely same-checksum.** Also
re-verified byte-identical to the `fold_ares` AF3 job in
[`analysis/afold/cofold_zn_three_species/SOURCE.md`](../../analysis/afold/cofold_zn_three_species/SOURCE.md).

```bash
# CRC64 C23878FFD0319265
curl "https://rest.uniprot.org/uniparc/search?query=checksum:C23878FFD0319265&format=json"
curl "https://rest.uniprot.org/uniparc/UPI00038F43E8?format=fasta"   # then diff
```

Two things follow. **Quote `WP_021575732` or `UPI00038F43E8`, never a UniProtKB
accession** — all six are dead links. And this explains why no *Escherichia* protein
carries a PF18870 annotation in UniProtKB (the query below returns zero rows): the
redundant *E. coli* entries were removed, not the protein. **UC5 aRES is not a strain
oddity — the identical 459-aa sequence is in ~150 *E. coli* assemblies.**

```
https://rest.uniprot.org/uniprotkb/search?query=(xref:pfam-PF18870)+AND+(organism_name:escherichia)&format=tsv
```

## `ares_kleb_pneumoniae.fasta` and `ares_enterobacter.fasta` — the two `docs/10` orthologs

**Identified and recorded 2026-10-14**, closing the README item that had them as VERIFY.
`docs/10` ran these as "Kleb_aRES" and "aRES_Echengduensis" without recording what they
were. Both sequences were recoverable from the AF3 job archives in
`analysis/afold/cofold_zn_three_species/` (the model coordinates carry the sequence), then
identified by the same CRC64 → UniParc → position-by-position comparison as above.

| file | `docs/10` label | aa | UniParc | RefSeq (active) | dead UniProtKB |
|---|---|---|---|---|---|
| `ares_kleb_pneumoniae.fasta` | Kleb_aRES | 443 | **UPI0006936E4D** | **WP_052450931** | `A0A169QRT4` |
| `ares_enterobacter.fasta` | aRES_Echengduensis | 459 | **UPI0003BED449** | **WP_023337974** | `V3RF86` |

- **`Kleb_aRES` is confirmed**: *Klebsiella pneumoniae*. The same sequence also occurs in
  *K. variicola* subsp. *variicola* (`XFQ07495`, `XGP45966`), so it is a
  *K. pneumoniae*-complex protein rather than one species'.
- ⚠️ **`aRES_Echengduensis`'s species is NOT confirmed, and the genus is *Enterobacter*,
  not *Escherichia*.** Every cross-reference is in the *Enterobacter cloacae* complex —
  *Enterobacter* sp. MGH 23 (`ESN20174`), *E. kobei* (`ELE9684022`, `EOG1790010`),
  *E. roggenkampii* (`MDU5045294`) — and RefSeq calls it "Enterobacter cloacae complex".
  So "E." almost certainly meant ***Enterobacter* chengduensis**, which is a member of
  that complex, and the label was right in genus after all. But *E. chengduensis* itself
  appears in none of the cross-references, and the identical sequence is shared across
  several complex members, so **the species cannot be recovered from the sequence** and
  the file is named for the genus. If the species matters, it has to come from whatever
  record the `docs/10` submission was built from.

These two are **not** in the `docs/11` panel. They are far closer to UC5 (71.9% and
44.5%) than anything available in UniProtKB with a PF18870 annotation, which is exactly
why they were worth identifying — and exactly why `docs/11` could not use them, having
fixed its selection filters on the annotated set before these were found.

How the sequences were recovered, for the record — `read_chain` returns the sequence
built from a model's own residues, which is why `homolog_compare.py` can take a structure
in place of a FASTA:

```bash
python3 - <<'EOF'
import os, sys, zipfile
sys.path.insert(0, "analysis/scripts")
from homolog_compare import read_chain
d = "analysis/afold/cofold_zn_three_species"
for z in ("fold_kleb_ares.zip", "fold_ares_echengduensis.zip"):
    with zipfile.ZipFile(os.path.join(d, z)) as zf:
        name = next(n for n in sorted(zf.namelist()) if n.endswith("model_0.cif"))
        open("/tmp/m.cif", "wb").write(zf.read(name))
    print(z, read_chain("/tmp/m.cif")[0])
EOF
```

---

## `pf18870_res_panel18.fasta` — the 18-ortholog test panel (`docs/11`)

**Retrieved 2026-10-14.** UniProtKB release **2026_03** (release date 02-September-2026,
read from the `X-UniProt-Release` response header). Domain coordinates from **InterPro**
(PF18870 = IPR041206, "HEPN/RES N-terminal domain 1"; PF08808 = RES domain), same date.

Candidate pool — proteins carrying **both** Pfam domains, which is what defines the
family this project is about (**213 hits**):

```
https://rest.uniprot.org/uniprotkb/search?query=(xref:pfam-PF18870)+AND+(xref:pfam-PF08808)&format=fasta&size=500
```

For reference, PF18870 alone returns 233 proteins at this release, so **213/233 (91%) of
all PF18870 proteins are RES-fused.** Domain coordinates, one call per protein:

```
https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/<ACCESSION>/
```

### Selection filters, fixed before any per-region identity was computed

The quantity under test in `docs/11` is per-region identity, so the panel had to be
chosen on criteria blind to it. In order:

| | filter | kept |
|---|---|---|
| F1 | PF18870 and PF08808 both matched, PF18870 N-terminal to PF08808 | 213/213 |
| F2 | not a UniProt sequence fragment | 213/213 |
| F3 | overall identity to UC5 ≥ 25% (the pool's 5th percentile is 24.9%; global alignment below ~25% is not trustworthy) | 199 |
| F4 | ≥ 60% of UC5 positions aligned in each of 1–48, 49–120, 126–447, so no region's identity rests on a handful of columns | 190 |
| F5 | one sequence per genus, the highest-identity member, so no genus can dominate | 126 |
| F6 | 18 taken at even rank intervals across the identity range, so the panel spans the family's divergence instead of clustering at the top | **18** |

F3 and F4 condition on *overall* divergence and on coverage, never on the identity of
one region against another. F4's effect is reported in `docs/11` §6 because coverage and
conservation are not independent: it drops 9 sequences and its one real risk — excluding
proteins that have deleted a region — is checked there with the indel-inclusive metric,
which needs no coverage filter.

### The panel

Identity is to UC5 over a global pairwise alignment (BLOSUM62, affine gaps). Pfam
columns are each protein's **own** numbering, as InterPro reports it.

| accession | organism | aa | id. to UC5 | PF18870 | PF08808 | phylum/clade |
|---|---|---|---|---|---|---|
| `A0A2M9WDD2` | *Pantoea rodasii* | 451 | 41.0% | 26–150 | 191–357 | Pseudomonadota |
| `B1Y0A5` | *Leptothrix cholodnii* SP-6 | 481 | 35.0% | 34–172 | 239–379 | Pseudomonadota |
| `A0A521DN47` | *Pedobacter westerhofensis* | 482 | 33.9% | 29–165 | 228–377 | FCB group |
| `A0A5N7MAN0` | *Microvirga tunisiensis* | 489 | 33.1% | 30–172 | 224–379 | Pseudomonadota |
| `R4WPB5` | *Caballeronia insecticola* | 499 | 32.3% | 144–268 | 307–472 | Pseudomonadota |
| `A0ABT1D458` | *Siccirubricoccus soli* | 394 | 31.1% | 36–163 | 201–365 | Pseudomonadota |
| `A0A2P7NRM1` | *Nitrosomonas supralitoralis* | 404 | 30.6% | 41–166 | 230–383 | Pseudomonadota |
| `A0ABV6FVC9` | *Fontibacter flavus* | 410 | 30.1% | 39–164 | 209–383 | FCB group |
| `A0A919Q5A9` | *Demequina activiva* | 498 | 29.5% | 127–245 | 315–455 | Actinomycetota |
| `A0A2S7X853` | *Aliivibrio sifiae* | 390 | 29.1% | 41–168 | 209–377 | Pseudomonadota |
| `A0ABY6A846` | *Thalassolituus hydrocarboniclasticus* | 450 | 28.6% | 41–166 | 232–374 | Pseudomonadota |
| `A0ABQ6Q009` | *Algoriphagus taiwanensis* | 419 | 28.3% | 38–168 | 234–380 | FCB group |
| `A0ABW1T0I8` | *Longivirga aurantiaca* | 393 | 27.9% | 40–163 | 202–360 | Actinomycetota |
| `A0ABU3NJY4` | *Halomonas saccharevitans* | 406 | 27.6% | 41–168 | 218–374 | Pseudomonadota |
| `A0A840RZ99` | *Inhella inkyongensis* | 398 | 27.0% | 40–165 | 222–376 | Pseudomonadota |
| `F0S5Z8` | *Pseudopedobacter saltans* DSM 12145 | 383 | 26.7% | 41–165 | 195–371 | FCB group |
| `A0A371JDE8` | *Lachnotalea glycerini* | 410 | 26.1% | 28–166 | 225–360 | Bacillota |
| `A0ABM7SVF0` | *Mycobacterium senriense* | 437 | 25.3% | 43–167 | 230–367 | Actinomycetota |

18 genera, 4 phyla/clades (Pseudomonadota 10, FCB 4, Actinomycetota 3, Bacillota 1).
Same table, machine-readable, with taxon IDs:
[`results/panel18_orthologs.tsv`](../../results/panel18_orthologs.tsv).

⚠️ **All 18 are UniProtKB *unreviewed* (TrEMBL) entries with protein existence
"Predicted"** — as are all 213 in the pool, and as is every PF18870 protein in the
database. None has experimental evidence of expression, let alone of anti-phage
function. They are gene predictions from genome assemblies. This is a property of the
whole family, not of the selection, but it bounds what the panel can support: it is
evidence about *sequence* evolution in a predicted gene family and nothing more.

## `pf18870_res_pool190.fasta` — the sensitivity set

The 190 survivors of F1–F4, i.e. the panel before genus de-duplication and rank
thinning. Same retrieval, same date. Used in `docs/11` §6 to show that the result does
not depend on the choice of 18, and to demonstrate the opposite failure: at n = 190 the
uncorrected test reports p < 0.0001 for a comparison the phylogenetic correction puts at
p = 0.09–0.30. Internal identity runs 19.5–100% — it contains identical sequences, which
is exactly why it is the sensitivity set and not the panel.

Accessions: the FASTA headers, and the `acc` column of
[`results/pool190_region_identity.tsv`](../../results/pool190_region_identity.tsv).

### Reproducing both files

```bash
# 213 dual-domain proteins (release is in the X-UniProt-Release response header)
curl -D headers.txt -o pool.fasta \
  "https://rest.uniprot.org/uniprotkb/search?query=(xref:pfam-PF18870)+AND+(xref:pfam-PF08808)&format=fasta&size=500"
```

Then apply F1–F6. The identity, coverage and taxonomy that F3–F6 need are the
`--identity-tsv` output of
[`analysis/scripts/homolog_compare.py`](../../analysis/scripts/homolog_compare.py) plus
the UniProt `lineage` field; the InterPro call above supplies F1.

A later UniProt release will not reproduce these 18 exactly — the pool grows and the even
rank spacing in F6 then lands elsewhere. The accessions above are the record; the filters
are the method.
