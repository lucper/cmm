SHELL       := /usr/bin/bash
.SHELLFLAGS := -o pipefail -c

CMM_TEST := ../../bin/cmm_test

.PHONY: usage
usage:
	@echo "Usage: make all"
	@echo "       INSTANCE=<instance name>"
	@echo "       FASTA=<fasta file>"
	@echo "       INTS=<interactions file>"
	@echo "       L=<motif length>"
	@echo "       D=<number of wildcards>"
	@echo "       K=<number of top motif pairs>"
	@echo "       THREADS=<number of threads>"
	@echo "       TRIAL=<trial number>"

.DELETE_ON_ERROR:

# One measurement (one line of cmm_test output) per file.
PERF := $(INSTANCE)/$(INSTANCE).l$(L)d$(D).k$(K).p$(THREADS).trial$(TRIAL).perf

.PHONY: all
all: $(PERF)

# The output is renamed only once cmm_test has finished, so that a killed run
# never leaves an incomplete measurement that would be skipped later.
$(PERF): $(FASTA) $(INTS)
	mkdir -p $(@D)
	$(CMM_TEST) $(FASTA) $(INTS) $(L) $(D) $(K) $(THREADS) /dev/null \
		> $@.tmp
	mv $@.tmp $@
