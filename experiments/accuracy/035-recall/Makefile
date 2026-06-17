SHELL := /usr/bin/bash

INPUT  := ../030-eval/$(DATASET)

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
	@echo "       C=<cutoff for recall in [0..100]>"

.DELETE_ON_ERROR:

.PHONY: all
all: $(DATASET)/$(DATASET).seq_slider.l$(L)d$(D).recall_h$(H)c$(C).trials$(R).min$(T).dat \
     $(DATASET)/$(DATASET).m_slider.l$(L)d$(D).recall_h$(H)c$(C).trials$(R).min$(T).dat

$(DATASET)/$(DATASET).%.l$(L)d$(d).recall_h$(H)c$(C).trials$(R).min$(T).dat: $(INPUT)/$(DATASET).%.l$(L)d$(D).agg.trials$(R).min$(T).k$(K_TOP).h$(H).dedup.tsv
	mkdir -p $(@D)
	LC_ALL=C awk -v t="$(C)" '{if (($$8 + 0) * 100 >= t) sum += 1; total += 1} END {if (total > 0) printf "%d %d %d%%\n", sum, total, (sum/total)*100; else print "0 0 NA"}' $< > $@
