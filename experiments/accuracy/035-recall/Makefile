SHELL := /usr/bin/bash

EVAL_INPUT  := ../030-eval/$(DATASET)
DEDUP_INPUT := ../025-dedup/$(DATASET)

.PHONY: usage
usage:
	@echo "Usage: make DATASET=<dataset name>"
	@echo "            R=<number of runs>"
	@echo "            T=<max convergence time>"
	@echo "            ELL=<motif length>"
	@echo "            D=<number of wildcards>"
	@echo "            K_TOP=<top motif pairs>"
	@echo "            <target>"

.DELETE_ON_ERROR:

$(DATASET)/$(DATASET).recall.l$(ELL)d$(D).trials$(R).min$(T).tsv: $(EVAL_INPUT)/$(DATASET).slider++.l$(ELL)d$(D).agg.trials$(R).min$(T).k$(K_TOP).h0.dedup.tsv $(DEDUP_INPUT)/$(DATASET).cmm.l$(ELL)d$(D).dedup.out
	mkdir -p $(@D)
	./compute_recall.sh $^ $@
