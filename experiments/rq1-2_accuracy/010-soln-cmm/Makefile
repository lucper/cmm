SHELL       := /usr/bin/bash
.SHELLFLAGS := -o pipefail -c

CMM_TEST := ../../bin/cmm_test

.PHONY: usage
usage:
	@echo "Usage: make all"
	@echo "       DATASET=<dataset name>"
	@echo "       FASTA=<fasta file>"
	@echo "       INTS=<interactions file>"
	@echo "       L=<motif length>"
	@echo "       D=<number of wildcards>"
	@echo "       K_TOP=<number of top motif pairs>"
	@echo "       THREADS=<number of threads>"

.DELETE_ON_ERROR:

OUT  := $(DATASET)/$(DATASET).cmm.l$(L)d$(D).k$(K_TOP).out
PERF := $(DATASET)/$(DATASET).cmm.l$(L)d$(D).k$(K_TOP).p$(THREADS).perf

.PHONY: all
all: $(PERF) $(OUT)

$(PERF): $(FASTA) $(INTS)
	mkdir -p $(@D)
	$(CMM_TEST) $(FASTA) $(INTS) $(L) $(D) $(K_TOP) $(THREADS) $(OUT).tmp \
		> $@
	if [ -f $(OUT) ]; then rm $(OUT).tmp; else mv $(OUT).tmp $(OUT); fi

$(OUT): $(FASTA) $(INTS)
	mkdir -p $(@D)
	$(CMM_TEST) $(FASTA) $(INTS) $(L) $(D) $(K_TOP) $(THREADS) $@ \
		> $(PERF)
