SHELL := /usr/bin/bash
.RECIPEPREFIX := >
.DELETE_ON_ERROR:

VENV_DIR    := venv
PYTHON      := $(VENV_DIR)/bin/python
PIP         := $(VENV_DIR)/bin/pip

# default parameters
PROGS       := cmm slider++
ED          := 05 10 15 20 25 30
L           := 8
D           := 3
K           := 10
S           := 100
DATASET     := yeast
TOP_K       := 50
OUT_PREFIX  := sampled_$(DATASET)_l$(L)_d$(D)_K$(K)_S$(S)
PROG_PREFIX := sampled_$(DATASET)_l$(L)_d$(D)_k$(TOP_K)_S$(S)
RESULTS     := $(foreach m,$(PROGS),results_$(m).txt)

usage:
> # DATASET ED SEQ INTS PROG
.PHONY: usage

clear:
> rm -rf *.png
.PHONY: clear

all: comparison_precision.png comparison_recall.png
.PHONY: all

comparison_precision.png comparison_recall.png &: $(RESULTS)
> cat $^ | $(PYTHON) plot.py

results_%.txt: $(foreach i,$(ED),$(PROG_PREFIX)_e$(i).%.out)
> rm -f $@
> for i in $(ED); do \
>   echo -n "$* 0.$${i} " >> $@; \
>   $(PYTHON) evaluate_precision_and_recall.py \
> 	  --predicted $(PROG_PREFIX)_e$${i}.$*.out \
>	  --k $(TOP_K) \
> 	  --truth $(OUT_PREFIX)_e$${i}.dat \
>	  --sequences $(OUT_PREFIX)_e$${i}.fa \
>	  --interactions $(OUT_PREFIX)_e$${i}.int >> $@; \
> done

$(OUT_PREFIX)_e%.fa $(OUT_PREFIX)_e%.int $(OUT_PREFIX)_e%.dat &: $(VENV_DIR)/bin/activate
> $(PYTHON) ppi_network_generator.py \
> 	--fasta $(SEQ) \
> 	--interactions $(INTS) \
> 	--l $(L) --d $(D) \
> 	--K $(K) --S $(S) --e 0.$* \
> 	--seed 42 \
> 	--out-fasta $(OUT_PREFIX)_e$*.fa \
> 	--out-interactions $(OUT_PREFIX)_e$*.int \
> 	--out-pairs $(OUT_PREFIX)_e$*.dat

$(VENV_DIR)/bin/activate: requirements.txt
> python3 -m venv $(VENV_DIR)
> $(PIP) install --upgrade pip
> $(PIP) install -r $<
