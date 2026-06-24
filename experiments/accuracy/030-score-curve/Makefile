SHELL       := /usr/bin/bash

INPUT       := ../025-dedup/$(DATASET)
EXACT_LABEL := Exact
EXACT_PROG  := cmm
METHODS_CSV := ../methods.csv
METHODS     := $(shell awk -F, '{print $$1}' $(METHODS_CSV))
WIDTH       := 3.8
HEIGHT      := 3.2

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

curve_legend.tex: ../methods.csv
	gnuplot -c plot_curve.gp legend $< curve_legend $(EXACT_LABEL) 7 1.1

$(DATASET)/curve.%.tex: $(DATASET)/%.dat
	gnuplot -c plot_curve.gp plot $< $(DATASET)/curve.$* $(EXACT_LABEL) $(WIDTH) $(HEIGHT)

$(DATASET)/$(DATASET).l$(L)d$(D).dedup.dat: ../methods.csv \
	$(INPUT)/$(DATASET).*.l$(L)d$(D).agg.trials$(R).min$(T).k$(K_TOP).dedup.out \
	$(INPUT)/$(DATASET).$(EXACT_PROG).l$(L)d$(D).k$(K_TOP).dedup.out
	mkdir -p $(@D)
	./glob.py $(EXACT_PROG) $^ > $@
