SHELL      := /usr/bin/bash

CMM_DEDUP  := ../cmm_dedup

## change these
DATA       := ../../string_data_cleaned
EXACT_DIR  := ../cmm_out

FASTA      := $(DATA)/$(DATASET).s700.cleaned.fa
INTS       := $(DATA)/$(DATASET).s700.cleaned.int

INPUT := ../020-agg/$(DATASET)

.PHONY: usage
usage:
	@echo "Usage: make DATASET=<dataset name>"
	@echo "            R=<number of runs>           (slider dedup only)"
	@echo "            T=<max convergence time>     (slider dedup only)"
	@echo "            ELL=<motif length>"
	@echo "            D=<number of wildcards>"
	@echo "            K_TOP=<top motif pairs>      (slider dedup only)"
	@echo "            <target>"

.DELETE_ON_ERROR:

## Dedup SLIDER's top-K aggregates.
$(DATASET)/$(DATASET).slider++.l$(ELL)d$(D).agg.trials$(R).min$(T).k$(K_TOP).dedup.out: $(INPUT)/$(DATASET).slider++.l$(ELL)d$(D).agg.trials$(R).min$(T).k$(K_TOP).out $(FASTA) $(INTS)
	mkdir -p $(@D)
	$(CMM_DEDUP) $(FASTA) $(INTS) $< 0 $@

## Dedup exact solution.
$(DATASET)/$(DATASET).cmm.l$(ELL)d$(D).k$(K_TOP).dedup.out: $(EXACT_DIR)/$(DATASET).cmm.l$(ELL)d$(D).k$(K_TOP).out $(FASTA) $(INTS)
	mkdir -p $(@D)
	$(CMM_DEDUP) $(FASTA) $(INTS) $< 0 $@
