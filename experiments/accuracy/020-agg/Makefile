SHELL := /usr/bin/bash

INPUT := ../010-soln/$(DATASET)

.PHONY: usage
usage:
	@echo "Usage: make DATASET=<dataset name>"
	@echo "            R=<number of runs>"
	@echo "            T=<max convergence time>"
	@echo "            ELL=<l>"
	@echo "            D=<d>"
	@echo "            K=<number of top motif pairs>"
	@echo "            K_TOP=<number of top motif pairs to keep>"
	@echo "            <target>"

.DELETE_ON_ERROR:

## RQ1 cells: R=5, T in {5, 10, 15, 20, 25, 30}
## RQ2 cells: T=5, R in {5, 10, 15, 20, 25, 30}
$(DATASET)/$(DATASET).slider++.l$(ELL)d$(D).agg.trials$(R).min$(T).k$(K_TOP).out: $(foreach r,$(shell seq 1 $(R)),$(INPUT)/$(DATASET).slider++.l$(ELL)d$(D).trial$(r).min$(T).k$(K).out)
	mkdir -p $(@D)
	cat $^ | sort -u | sort -k3 -gr | awk -v n=$(K_TOP) 'NR<=n' > $@
