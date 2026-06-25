SHELL       := /usr/bin/bash

INPUT       := ../030-eval/$(DATASET)
EXACT_LABEL := Exact
EXACT_PROG  := cmm
METHODS_CSV := ../methods.csv
METHODS     := $(shell awk -F, '{print $$1}' $(METHODS_CSV))

## dimensions for plot and legend
PLOT_W      := 3.8
PLOT_H      := 3.2
LEGEND_W    := 7
LEGEND_H    := 1.1

THR         := 0.7

.PHONY: usage
usage:
	@echo "Usage: make all"
	@echo "       DATASET=<dataset name>"
	@echo "       R=<number of runs>"
	@echo "       T=<max convergence time>"
	@echo "       L=<motif length>"
	@echo "       D=<number of wildcards>"
	@echo "       K_TOP=<top motif pairs>"
	@echo "       THR=<threshold line for plot>"

.DELETE_ON_ERROR:

.PHONY: all
all: hist.$(DATASET).l$(L)d$(D).dedup.tex hist_legend.tex

hist_legend.tex: $(METHODS_CSV)
	gnuplot -c plot_hist.gp legend \
		$< \
		hist_legend \
		$(EXACT_PROG) \
		$(LEGEND_W) $(LEGEND_H)

hist.%.l$(L)d$(D).dedup.tex: $(INPUT)/%.*.l$(L)d$(D).agg.trials$(R).min$(T).k$(K_TOP).h$(H).dedup.out
	gnuplot -c plot_hist.gp plot \
		$(METHODS_CSV) \
		hist.$*.l$(L)d$(D).dedup \
		$(EXACT_PROG) \
		$(PLOT_W) $(PLOT_H) \
		$^ \
		$(THR)
