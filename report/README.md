# Findings report

`aRES_findings.tex` — summary of the structural findings so far. Build with
`make`, or `pdflatex aRES_findings.tex` twice.

Needs: `pgfplots`, `multirow`, `siunitx`, `mdframed`, `booktabs`, `hyperref`.

## Figures

| File | Status | Source |
|---|---|---|
| `figures/fig_domain_map.tex` | generated (TikZ) | boundaries from `docs/05` §1.1, `docs/07` §1 |
| `figures/fig_plddt.tex` + `plddt.dat` | generated (pgfplots) | B-factor column of the monomer model |
| `figures/fig_cartoon.png` | **supply** | PyMOL, `analysis/scripts/ares_domain_overview.pml` scene `F1_domains` |
| `figures/fig_zincsite.png` | **supply** | PyMOL scene `F4_zincsite`, or the Coot view |

The two TikZ figures are vector and rebuild from source, so they track the data.
Regenerate `plddt.dat` from a different model with:

```bash
python3 - <<'PY'
import statistics
per = {}
for ln in open("model.pdb"):
    if ln.startswith("ATOM  "):
        per.setdefault(int(ln[22:26]), []).append(float(ln[60:66]))
with open("figures/plddt.dat", "w") as fh:
    fh.write("resi plddt\n")
    for r in sorted(per):
        fh.write(f"{r} {statistics.mean(per[r]):.2f}\n")
PY
```
