SHELL       := /usr/bin/bash

INPUT       := ../025-dedup/$(DATASET)
EXACT_PROG  := cmm
METHODS_CSV := ../methods.csv
METHODS     := $(shell awk -F, '{print $$1}' $(METHODS_CSV))

## dimensions of the plot area (the same in all plots of the experiments) and of the legend, in cm
PLOT_W      := 5.0
PLOT_H      := 3.6
LEGEND_W    := 7
LEGEND_H    := 1.0

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
all: curve.$(DATASET).l$(L)d$(D).dedup.tex curve_legend.tex

curve_legend.tex:
	gnuplot -c plot_curve.gp legend \
		$(METHODS_CSV) \
		curve_legend \
		$(EXACT_PROG) \
		$(LEGEND_W) $(LEGEND_H)

curve.%.tex: %.dat
	gnuplot -c plot_curve.gp plot \
		$(METHODS_CSV) \
		$< \
		curve.$* \
		$(EXACT_PROG) \
		$(PLOT_W) $(PLOT_H)

$(DATASET).l$(L)d$(D).dedup.dat: \
	$(INPUT)/$(DATASET).*.l$(L)d$(D).agg.trials$(R).min$(T).k$(K_TOP).dedup.out \
	$(INPUT)/$(DATASET).$(EXACT_PROG).l$(L)d$(D).k$(K_TOP).dedup.out
	./glob.py $(EXACT_PROG) $(METHODS_CSV) $^ > $@
