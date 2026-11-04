SHELL       := /usr/bin/bash
.SHELLFLAGS := -o pipefail -c

CMM_TEST_MPI := ../../bin/cmm_test_mpi
CMM_TEST     := ../../bin/cmm_test

.PHONY: usage
usage:
	@echo "Usage: make -k all, inside a Slurm job on NODES whole nodes (without -j: a run takes them all)"
	@echo "       or make -k no-mpi, to run cmm_test instead, on 1 whole node, with THREADS threads"
	@echo "       DATASET=<dataset name>"
	@echo "       FASTA=<fasta file>"
	@echo "       INTS=<interactions file>"
	@echo "       L=<motif length>"
	@echo "       D=<number of wildcards>"
	@echo "       F=<support function>"
	@echo "       K=<number of top motif pairs>"
	@echo "       NODES=<number of nodes>"
	@echo "       LAYOUTS=<list of processes per node>x<threads per process>"
	@echo "       THREADS=<number of threads, for no-mpi>"
	@echo "       TRIALS=<number of trials>"

.DELETE_ON_ERROR:

PREFIX := $(DATASET)/$(DATASET).layouts.l$(L)d$(D).$(F).k$(K).n$(NODES)

# One file per layout and trial, in turn. Those that exist are skipped.
PERFS := $(foreach trial,$(shell seq 1 $(TRIALS)),$(foreach layout,$(LAYOUTS),$(PREFIX).$(layout).trial$(trial).perf))

# The current code without MPI, as if it were one more layout, so that 020-concat treats it as the others.
NO_MPI_PERFS := $(foreach trial,$(shell seq 1 $(TRIALS)),$(PREFIX).nompi$(THREADS).trial$(trial).perf)

.PHONY: all no-mpi
all: $(PERFS)
no-mpi: $(NO_MPI_PERFS)

# In a recipe, $* is <layout>.trial<trial>.
layout    = $(word 1,$(subst ., ,$*))
per_node  = $(word 1,$(subst x, ,$(layout)))
threads   = $(word 2,$(subst x, ,$(layout)))
processes = $(shell echo $$(( $(per_node) * $(NODES) )))

# One measurement (one line of cmm_test_mpi output) per file. The outputs are renamed only once
# cmm_test_mpi has finished, so that a job killed at its time limit never leaves an incomplete
# measurement that would be skipped later. block:block keeps the cores of a process together, where
# Slurm's default may spread them over the sockets; verbose writes them to the log of the job.
$(PERFS): $(PREFIX).%.perf: $(FASTA) $(INTS)
	mkdir -p $(@D)
	srun --nodes=$(NODES) --ntasks=$(processes) --ntasks-per-node=$(per_node) --cpus-per-task=$(threads) \
	     --distribution=block:block --cpu-bind=verbose \
	     $(CMM_TEST_MPI) $(FASTA) $(INTS) $(L) $(D) $(K) $(threads) $(@:.perf=.out).tmp $(F) > $@.tmp 2> $(@:.perf=.stderr)
	grep '^p[0-9]' $(@:.perf=.stderr) > $(@:.perf=.tsv)
	mv $(@:.perf=.out).tmp $(@:.perf=.out)
	mv $@.tmp $@

$(NO_MPI_PERFS): $(FASTA) $(INTS)
	mkdir -p $(@D)
	$(CMM_TEST) $(FASTA) $(INTS) $(L) $(D) $(K) $(THREADS) $(@:.perf=.out).tmp $(F) > $@.tmp
	mv $(@:.perf=.out).tmp $(@:.perf=.out)
	mv $@.tmp $@
