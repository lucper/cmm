SHELL     := /usr/bin/bash

CMM_EVAL  := ../cmm_eval

## change this
DATA      := ../../string_data_cleaned
FASTA     := $(DATA)/$(DATASET).s700.cleaned.fa
INTS      := $(DATA)/$(DATASET).s700.cleaned.int

INPUT     := ../020-agg/$(DATASET)

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
$(DATASET)/$(DATASET).slider++.l$(ELL)d$(D).agg.trials$(R).min$(T).k$(K_TOP).h$(H).tsv: $(INPUT)/$(DATASET).slider++.l$(ELL)d$(D).agg.trials$(R).min$(T).k$(K_TOP).out
	mkdir -p $(@D)
	$(CMM_EVAL) $(FASTA) $(INTS) $< $(EXACT_SOL) $(H) $@ /dev/null
