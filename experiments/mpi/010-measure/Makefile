SHELL       := /usr/bin/bash
.SHELLFLAGS := -o pipefail -c

CMM_TEST_MPI := ../../bin/cmm_test_mpi

.PHONY: usage
usage:
	@echo "Usage: make all, inside a Slurm job that has the cores of LAYOUT"
	@echo "       DATASET=<dataset name>"
	@echo "       FASTA=<fasta file>"
	@echo "       INTS=<interactions file>"
	@echo "       L=<motif length>"
	@echo "       D=<number of wildcards>"
	@echo "       K=<number of top motif pairs>"
	@echo "       LAYOUT=<processes>x<threads per process>"
	@echo "       TRIAL=<trial number>"

.DELETE_ON_ERROR:

NAME := $(DATASET)/$(DATASET).layouts.l$(L)d$(D).k$(K).$(LAYOUT).trial$(TRIAL)
PERF := $(NAME).perf
OUT  := $(NAME).out

PROCESSES := $(word 1,$(subst x, ,$(LAYOUT)))
THREADS   := $(word 2,$(subst x, ,$(LAYOUT)))

.PHONY: all
all: $(PERF)

# One measurement (one line of cmm_test_mpi output) per file. The outputs are renamed only once
# cmm_test_mpi has finished, so that a job killed at its time limit never leaves an incomplete
# measurement that would be skipped later. block:block keeps the cores of a process together, where
# Slurm's default may spread them over the sockets.
$(PERF): $(FASTA) $(INTS)
	mkdir -p $(@D)
	srun --nodes=1 --ntasks=$(PROCESSES) --ntasks-per-node=$(PROCESSES) --cpus-per-task=$(THREADS) \
	     --distribution=block:block \
	     $(CMM_TEST_MPI) $(FASTA) $(INTS) $(L) $(D) $(K) $(THREADS) $(OUT).tmp x2 > $@.tmp
	mv $(OUT).tmp $(OUT)
	mv $@.tmp $@
