SHELL       := /usr/bin/bash
.SHELLFLAGS := -o pipefail -c

INPUT      := ../020-concat
PANELS     := ./make_panels.py
PLOT       := plot_efficiency.gp
ANNOTATE   := 8:4 8:7 5:2 5:4   # (ell:d) pairs whose time points are labelled with max_rank

## dimensions of the plot area (the same in all plots of the experiments) and of the legend, in cm
PLOT_W     := 5.0
PLOT_H     := 3.6
LEGEND_W   := 13
LEGEND_H   := 1.0

.PHONY: usage
usage:
	@echo "Usage: make all"
	@echo "       TSV=<TSV name in 020-concat (rq3.*, rq4.*, or rq5.*)>"

.DELETE_ON_ERROR:

.PHONY: all
all: $(TSV:.tsv=.plotted)

# make_panels.py splits the TSV into one file per panel, <TSV name without .tsv>_*.dat, and lists
# them in a .panels file; each panel file gets a time plot, a memory plot, and a legend, as .tex and
# .pdf files of the same name with _time, _memory, and _legend. The .plotted file marks that they
# are up to date.
%.plotted: $(INPUT)/%.tsv $(PANELS) $(PLOT)
	set -e; rq=$*; rq=$${rq%%.*}; \
	$(PANELS) $$rq $< $* "$(ANNOTATE)" > $*.panels; \
	for dat in $$(cat $*.panels); do \
		for mode in time memory; do \
			gnuplot -c $(PLOT) $$mode $$rq $$dat $${dat%.dat}_$$mode $(PLOT_W) $(PLOT_H); \
		done; \
		gnuplot -c $(PLOT) legend $$rq $$dat $${dat%.dat}_legend $(LEGEND_W) $(LEGEND_H); \
	done
	touch $@
