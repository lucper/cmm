SHELL       := /usr/bin/bash

INPUT       := ../025-dedup/$(DATASET)
EXACT_PROG  := cmm
METHODS_CSV := ../methods.csv
METHODS     := $(shell awk -F, '{print $$1}' $(METHODS_CSV))

## dimensions for plot and legend
PLOT_W      := 3.8
PLOT_H      := 3.2
LEGEND_W    := 7
LEGEND_H    := 1.1

.PHONY: usage
usage:
	@echo "Usage: make all"
	@echo "       DATASET=<dataset name>"
	@echo "       R=<number of runs>"
	@echo "       T=<max convergence time>"
	@echo "       L=<motif length>"
	@echo "       D=<number of wildcards>"
	@echo "       K_TOP=<top motif pairs>"

.DELETE_ON_ERROR:

.PHONY: all
all: $(DATASET)/curve.$(DATASET).l$(L)d$(D).dedup.tex curve_legend.tex

curve_legend.tex: $(METHODS_CSV)
	gnuplot -c plot_curve.gp legend \
		$(METHODS_CSV) \
		curve_legend \
		$(EXACT_PROG) \
		$(LEGEND_W) $(LEGEND_H)

$(DATASET)/curve.%.tex: $(DATASET)/%.dat
	gnuplot -c plot_curve.gp plot \
		$(METHODS_CSV) \
		$< \
		$(DATASET)/curve.$* \
		$(EXACT_PROG) \
		$(PLOT_W) $(PLOT_H)

$(DATASET)/$(DATASET).l$(L)d$(D).dedup.dat: \
	$(INPUT)/$(DATASET).*.l$(L)d$(D).agg.trials$(R).min$(T).k$(K_TOP).dedup.out \
	$(INPUT)/$(DATASET).$(EXACT_PROG).l$(L)d$(D).k$(K_TOP).dedup.out
	mkdir -p $(@D)
	./glob.py $(EXACT_PROG) $(METHODS_CSV) $^ > $@
