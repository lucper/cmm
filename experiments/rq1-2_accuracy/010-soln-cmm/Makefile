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

.PHONY: all
all: $(DATASET)/$(DATASET).cmm.l$(L)d$(D).k$(K_TOP).out

$(DATASET)/$(DATASET).cmm.l$(L)d$(D).k$(K_TOP).out: $(FASTA) $(INTS)
	mkdir -p $(@D)
	$(CMM_TEST) $(FASTA) $(INTS) $(L) $(D) $(K_TOP) $(THREADS) $@ \
		> $(@:.out=.perf)
