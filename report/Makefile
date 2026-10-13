# Build the findings report. Needs a TeX distribution with pgfplots and multirow
# (TeX Live full, MiKTeX, or MacTeX). Nothing here is generated in the repo.
TEX  = pdflatex
MAIN = aRES_findings

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex figures/fig_domain_map.tex figures/fig_plddt.tex figures/plddt.dat
	$(TEX) -interaction=nonstopmode $(MAIN).tex
	$(TEX) -interaction=nonstopmode $(MAIN).tex

clean:
	rm -f *.aux *.log *.out *.toc

distclean: clean
	rm -f $(MAIN).pdf

.PHONY: all clean distclean
