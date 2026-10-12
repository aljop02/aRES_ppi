# MD Plan: T7 Polymerase D106G / S122P and aRES Recognition

Oct 10, 2026 · @Felix

## Question and hypotheses

The simulations decide whether S122P and D106G act through the conformation of one loop in the T7 polymerase, through a direct contact to aRES, or through neither. They cannot test activation itself; they test what each mutation does to local structure, dynamics and contact persistence.

| Hypothesis | Mechanism | Predicted MD signature |
| --- | --- | --- |
| H1: loop conformation | Asp106 carboxylate caps the Gly121 N-H. Pro122 has no N-H and restricts phi(122) and psi(121). D106G removes the cap and adds backbone flexibility. | Both mutants show lower Asp106-Gly121 H-bond occupancy, shifted phi/psi at 121/122 or higher loop RMSF than WT. S122A stays WT-like. |
| H2: direct contact | In the cofold, Ser122 OG donates to Asp104 of aRES (3.06 Å). Pro cannot. | Loop stays WT-like in the S122P monomer. In the complex, Ser122-Asp104 occupancy is high in WT and zero in S122P, and S122A matches S122P. |
| H3: neither | The mutations destabilize the fold or act at another site (thioredoxin binding, or the cofold interface is simply wrong). | Local metrics unchanged, or changes far from the loop. The WT complex interface is itself unstable. Move to the wet-lab tests. |

A null result is still informative: if WT and mutants look alike and the WT interface is unstable, the cofold interface is not discriminating and the wet-lab tests in the last section carry the weight.

## Systems to simulate

Simulate the polymerase as it exists in the infected cell, gp5 bound to E. coli thioredoxin, and build every starting model from AlphaFold coordinates, because the crystal structure does not contain the loop. [1T7P](https://rcsb.org/structure/1t7p), the 2.2 Å structure of gp5 with thioredoxin, primer-template and dNTP, [carries a deletion of residues 118-123](https://ibdc.dbtindia.gov.in/isda/summary/1T7P). Gly121 and Ser122 are therefore absent from it, and the Asp106-Gly121 geometry exists only in your AF models. Use 1T7P to check the rest of the fold, not the loop. The surrounding stretch, residues 120-150, is [described as unique to T7 polymerase](https://richardson.med.harvard.edu/publications_190-231/205.pdf) when compared with the Klenow fragment, so expect it to be flexible and poorly templated.

| Tier | Systems | Contents | Purpose |
| --- | --- | --- | --- |
| 1 | gp5 WT, D106G, S122P | gp5 + thioredoxin, no DNA, mutants made in silico from the same AF models | Tests H1: loop conformation and the Asp106-Gly121 bond with no partner present |
| 2 | aRES + gp5 WT, D106G, S122P | The cofold complex as predicted, 1:1 | Tests H2: persistence of the Ser122-Asp104 contact and whether the WT interface is stable at all |
| 3 | gp5 D106A, D106N, S122A | Same as Tier 1 | Separates side-chain loss from flexibility (106) and from backbone restriction (122). Run if Tier 1 is ambiguous or supports H1 |
| 4 (optional) | gp5 from E. coli UC5 (activator) and E. chengduensis (non-activator) | AF models, same as Tier 1 | Checks whether the natural non-activator shows the same loop signature as the mutants |

Setup, before any production run:

1. Generate five AF models (different seeds) of gp5 + thioredoxin and of the aRES complex. Keep models whose 100-125 loop differs and record pLDDT for residues 100-125.
2. Superimpose each model on 1T7P outside residues 118-123 and record the Cα RMSD. If it is high, find out why before going on.
3. Make each mutant by in silico substitution of the same WT model (PyMOL or Rosetta), then minimize. Mutants then start identical to WT and differ only by the substitution.
4. Assign protonation states with PROPKA at the assay pH. Check Asp104, Asp106 and any His within 10 Å of the loop.
5. Confirm the chain assignment in the cofold files. This plan assumes residues 106, 121 and 122 belong to gp5 and Asp104 to aRES.
