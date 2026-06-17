SHELL       := /usr/bin/bash
.SHELLFLAGS := -o pipefail -c

CMM_DEDUP   := ../cmm_dedup

## change
EXACT_DIR   := ../cmm_out

INPUT       := ../020-agg/$(DATASET)

.PHONY: usage
usage:
	@echo "Usage: make [slider|cmm]"
	@echo "       DATASET=<dataset name>"
	@echo "       FASTA=<fasta file>"
	@echo "       INTS=<interactions file>"
	@echo "       R=<number of runs>           (SLIDER dedup only)"
	@echo "       T=<max convergence time>     (SLIDER dedup only)"
	@echo "       L=<motif length>"
	@echo "       D=<number of wildcards>"
	@echo "       K_TOP=<top motif pairs>      (SLIDER dedup only)"

.DELETE_ON_ERROR:

.PHONY: slider
slider: $(DATASET)/$(DATASET).seq_slider.l$(L)d$(D).agg.trials$(R).min$(T).k$(K_TOP).dedup.out \
	    $(DATASET)/$(DATASET).m_slider.l$(L)d$(D).agg.trials$(R).min$(T).k$(K_TOP).dedup.out

.PHONY: cmm
cmm: $(DATASET)/$(DATASET).cmm.l$(L)d$(D).k$(K_TOP).dedup.out

## Dedup SLIDER's top-K aggregates.
$(DATASET)/$(DATASET).%.l$(L)d$(D).agg.trials$(R).min$(T).k$(K_TOP).dedup.out: $(FASTA) $(INTS) $(INPUT)/$(DATASET).%.l$(L)d$(D).agg.trials$(R).min$(T).k$(K_TOP).out
	mkdir -p $(@D)
	$(CMM_DEDUP) $^ 0 /dev/stdout | LC_ALL=C sort -k3 -gr > $@

## Dedup exact solution.
$(DATASET)/$(DATASET).cmm.l$(L)d$(D).k$(K_TOP).dedup.out: $(FASTA) $(INTS) $(EXACT_DIR)/$(DATASET).cmm.l$(L)d$(D).k$(K_TOP).out
	mkdir -p $(@D)
	$(CMM_DEDUP) $^ 0 /dev/stdout | LC_ALL=C sort -k3 -gr > $@
