SHELL := /usr/bin/bash

INPUT := ../025-dedup/$(DATASET)

.PHONY: usage
usage:
	@echo "Usage: make all"
	@echo "       DATASET=<dataset name>"
	@echo "       R=<number of runs>"
	@echo "       T=<max convergence time>"
	@echo "       L=<l>"
	@echo "       D=<d>"
	@echo "       K_TOP=<number of top motif pairs to keep>"

.DELETE_ON_ERROR:

.PHONY: all
all: $(DATASET)/$(DATASET).seq_slider.l$(L)d$(D).normscore.trials$(R).min$(T).k$(K_TOP).dedup.dat \
	 $(DATASET)/$(DATASET).m_slider.l$(L)d$(D).normscore.trials$(R).min$(T).k$(K_TOP).dedup.dat

## The awk command assumes the files are sorted from highest to lowest scoring motifs.
## It then filters the motif pairs that do not have a corresponding one in the other file, effectively truncating the size to min(|S_A|,|S_B|).
$(DATASET)/$(DATASET).%.l$(L)d$(D).normscore.trials$(R).min$(T).k$(K_TOP).dedup.dat: $(INPUT)/$(DATASET).%.l$(L)d$(D).agg.trials$(R).min$(T).k$(K_TOP).dedup.out $(INPUT)/$(DATASET).cmm.l$(L)d$(D).k$(K_TOP).dedup.out
	mkdir -p $(@D)
	LC_ALL=C paste $^ -d' ' | awk 'NF==6 && $$6+0>0 {print $$3, $$6, ($$3/$$6) * 100}' > $@
