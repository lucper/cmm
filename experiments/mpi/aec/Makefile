SHELL       := /usr/bin/bash
.SHELLFLAGS := -o pipefail -c

# A separate check of the experiment of the layouts: the code of the tag ALENEX27_AEC, on one node, to compare its
# running time and memory with the MPI version. Nothing else depends on this directory.

TAG      := ALENEX27_AEC
CMM_TEST := source/experiments/bin/cmm_test

.PHONY: usage
usage:
	@echo "Usage: make build"
	@echo "       make -k run, inside a Slurm job on 1 whole node"
	@echo "       make merge"
	@echo "       DATASET=<dataset name>"
	@echo "       FASTA=<fasta file>, for run"
	@echo "       INTS=<interactions file>, for run"
	@echo "       L=<motif length>"
	@echo "       D=<number of wildcards>"
	@echo "       F=<support function>; the code of the tag only has x2"
	@echo "       K=<number of top motif pairs>"
	@echo "       THREADS=<number of threads>, for run"
	@echo "       TRIALS=<number of trials>, for run"

.DELETE_ON_ERROR:

NAME  := $(DATASET).layouts.l$(L)d$(D).$(F).k$(K)
PERFS := $(foreach trial,$(shell seq 1 $(TRIALS)),$(DATASET)/$(NAME).n1.aec$(THREADS).trial$(trial).perf)

.PHONY: build run merge
build: $(CMM_TEST)
run: $(PERFS)
merge: $(NAME).tsv

# The code of the tag has its own layout and programs, so it is kept apart. Once built, it stays: to build it
# again, on another partition, remove source/.
$(CMM_TEST):
	mkdir -p source
	git -C ../../.. archive $(TAG) | tar -x -C source
	$(MAKE) -C source/experiments/src ../bin/cmm_test

# One measurement per file, renamed only once the run has finished, as in ../010-measure.
$(PERFS): $(FASTA) $(INTS)
	[ "$(F)" = x2 ]
	mkdir -p $(@D)
	$(CMM_TEST) $(FASTA) $(INTS) $(L) $(D) $(K) $(THREADS) $(@:.perf=.out).tmp > $@.tmp
	mv $(@:.perf=.out).tmp $(@:.perf=.out)
	mv $@.tmp $@

# What ../020-concat/$(NAME).tsv has, if it is there yet, then the runs of the tag in the same format: the nodes, the
# layout and the trial, then the columns of cmm_test. Column 10 of the tag, the prunings, has no match, so it goes.
.PHONY: FORCE
$(NAME).tsv: FORCE
	cat $(wildcard ../020-concat/$(NAME).tsv) /dev/null > $@.tmp
	for perf in $(wildcard $(DATASET)/$(NAME).n1.aec*.perf); do \
	    awk -v OFS='\t' 'NF != 14 {print "FAILED: " FILENAME " has " NF " columns, not 14"; exit 1} \
	                     {n = split(FILENAME, parts, "."); trial = parts[n - 1]; sub(/^trial/, "", trial); \
	                      line = "1" OFS parts[n - 2] OFS trial; \
	                      for (i = 1; i <= NF; i++) if (i != 10) line = line OFS $$i; print line}' $$perf > $@.line \
	        || { cat $@.line; rm $@.line $@.tmp; exit 1; }; \
	    cat $@.line >> $@.tmp; rm $@.line; \
	done
	mv $@.tmp $@
