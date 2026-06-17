SHELL    := /usr/bin/bash

CMM_EVAL := ../cmm_eval

## change this
DATA     := ../../string_data_cleaned
FASTA    := $(DATA)/$(DATASET).s700.cleaned.fa
INTS     := $(DATA)/$(DATASET).s700.cleaned.int

INPUT    := ../025-dedup/$(DATASET)

.PHONY: usage
usage:
	@echo "Usage: make all"
	@echo "       DATASET=<dataset name>"
	@echo "       R=<number of runs>"
	@echo "       T=<max convergence time>"
	@echo "       L=<l>"
	@echo "       D=<d>"
	@echo "       K_TOP=<number of top motif pairs to keep>"
	@echo "       H=<similarity parameter>"

.DELETE_ON_ERROR:

.PHONY: all
all: $(DATASET)/$(DATASET).seq_slider.l$(L)d$(D).agg.trials$(R).min$(T).k$(K_TOP).h$(H).dedup.tsv \
	 $(DATASET)/$(DATASET).m_slider.l$(L)d$(D).agg.trials$(R).min$(T).k$(K_TOP).h$(H).dedup.tsv

$(DATASET)/$(DATASET).%.l$(L)d$(D).agg.trials$(R).min$(T).k$(K_TOP).h$(H).dedup.tsv: $(FASTA) $(INTS) $(INPUT)/$(DATASET).%.l$(L)d$(D).agg.trials$(R).min$(T).k$(K_TOP).dedup.out $(INPUT)/$(DATASET).cmm.l$(L)d$(D).k$(K_TOP).dedup.out
	mkdir -p $(@D)
	$(CMM_EVAL) $^ $(H) /dev/null $@
