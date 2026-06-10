SHELL       := /usr/bin/bash

CMM_EVAL    := ../cmm_eval

## change this
DATA        := ../../string_data_cleaned
FASTA       := $(DATA)/$(DATASET).s700.cleaned.fa
INTS        := $(DATA)/$(DATASET).s700.cleaned.int

AGG_INPUT   := ../020-agg/$(DATASET)
DEDUP_INPUT := ../025-dedup/$(DATASET)

.PHONY: usage
usage:
	@echo "Usage: make DATASET=<dataset name>"
	@echo "            R=<number of runs>"
	@echo "            T=<max convergence time>"
	@echo "            ELL=<l>"
	@echo "            D=<d>"
	@echo "            K_TOP=<number of top motif pairs to keep>"
	@echo "            H=<similarity parameter>"
	@echo "            EXACT_SOL=<path to exact solution file>"
	@echo "            <target>"

.DELETE_ON_ERROR:

$(DATASET)/$(DATASET).slider++.l$(ELL)d$(D).agg.trials$(R).min$(T).k$(K_TOP).h$(H).tsv: $(AGG_INPUT)/$(DATASET).slider++.l$(ELL)d$(D).agg.trials$(R).min$(T).k$(K_TOP).out $(EXACT_SOL)
	mkdir -p $(@D)
	$(CMM_EVAL) $(FASTA) $(INTS) $^ $(H) $@ /dev/null

$(DATASET)/$(DATASET).slider++.l$(ELL)d$(D).agg.trials$(R).min$(T).k$(K_TOP).h$(H).dedup.tsv: $(DEDUP_INPUT)/$(DATASET).slider++.l$(ELL)d$(D).agg.trials$(R).min$(T).k$(K_TOP).dedup.out $(DEDUP_INPUT)/$(DATASET).cmm.l$(ELL)d$(D).k$(K_TOP).dedup.out
	mkdir -p $(@D)
	$(CMM_EVAL) $(FASTA) $(INTS) $^ $(H) $@ /dev/null
