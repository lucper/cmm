SHELL       := /usr/bin/bash
.SHELLFLAGS := -o pipefail -c

INPUT       := ../010-soln/$(DATASET)

.PHONY: usage
usage:
	@echo "Usage: make all"
	@echo "       DATASET=<dataset name>"
	@echo "       R=<number of runs>"
	@echo "       T=<max convergence time>"
	@echo "       L=<l>"
	@echo "       D=<d>"
	@echo "       K=<number of top motif pairs>"
	@echo "       K_TOP=<number of top motif pairs to keep>"

.DELETE_ON_ERROR:

.PHONY: all
all: $(DATASET)/$(DATASET).seq_slider.l$(L)d$(D).agg.trials$(R).min$(T).k$(K_TOP).out \
     $(DATASET)/$(DATASET).m_slider.l$(L)d$(D).agg.trials$(R).min$(T).k$(K_TOP).out

## RQ1 cells: R=5, T in {5, 10, 15, 20, 25, 30}
## RQ2 cells: T=5, R in {5, 10, 15, 20, 25, 30}
$(DATASET)/$(DATASET).%.l$(L)d$(D).agg.trials$(R).min$(T).k$(K_TOP).out: $(foreach r,$(shell seq 1 $(R)),$(INPUT)/$(DATASET).%.l$(L)d$(D).trial$(r).min$(T).k$(K).out)
	mkdir -p $(@D)
	export LC_ALL=C; cat $^ | sort -u | sort -k3 -gr | awk -v n=$(K_TOP) 'NR<=n' > $@
