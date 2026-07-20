SHELL       := /usr/bin/bash

INPUT       := ../030-eval/$(DATASET)
EXACT_LABEL := Exact
EXACT_PROG  := cmm
METHODS_CSV := ../methods.csv
METHODS     := $(shell awk -F, '{print $$1}' $(METHODS_CSV))

## dimensions for plot and legend
PLOT_W      := 3.5
PLOT_H      := 3.0
LEGEND_W    := 7
LEGEND_H    := 1.0

THR         := 0 # default to no dashed line

.PHONY: usage
usage:
	@echo "Usage: make all"
	@echo "       DATASET=<dataset name>"
	@echo "       R=<number of runs>"
	@echo "       T=<max convergence time>"
	@echo "       L=<motif length>"
	@echo "       D=<number of wildcards>"
	@echo "       K_TOP=<top motif pairs>"
	@echo "       H=<similarity parameter>"
	@echo "       THR=<threshold line for plot>"

.DELETE_ON_ERROR:

.PHONY: all
all: density.$(DATASET).l$(L)d$(D).h$(H).dedup.tex density_legend.tex

density_legend.tex:
	gnuplot -c plot_density.gp legend \
		$(METHODS_CSV) \
		density_legend \
		$(EXACT_PROG) \
		$(LEGEND_W) $(LEGEND_H)

density.%.l$(L)d$(D).h$(H).dedup.tex: $(INPUT)/%.*.l$(L)d$(D).agg.trials$(R).min$(T).k$(K_TOP).h$(H).dedup.tsv
	gnuplot -c plot_density.gp plot \
		$(METHODS_CSV) \
		density.$*.l$(L)d$(D).h$(H).dedup \
		$(EXACT_PROG) \
		$(PLOT_W) $(PLOT_H) \
		$^ \
		$(THR)
