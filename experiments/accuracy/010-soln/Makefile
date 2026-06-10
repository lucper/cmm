SHELL   := /usr/bin/bash

SLIDER  := java -cp ../000-slider/SliderLight/dist/SliderLight.jar Framework.Framework

## change this
DATA    := ../../string_data_cleaned
FASTA   := $(DATA)/$(DATASET).s700.cleaned.fa
INTS    := $(DATA)/$(DATASET).s700.cleaned.int

.PHONY: usage
usage:
	@echo "Usage: make DATASET=<dataset name>"
	@echo "            I=<run index>"
	@echo "            T=<max convergence time>"
	@echo "            ELL=<motif length>"
	@echo "            D=<number of wildcards>"
	@echo "            K=<number of top motif pairs>"
	@echo "            <target>"

.DELETE_ON_ERROR:

$(DATASET)/$(DATASET).slider++.l$(ELL)d$(D).trial$(I).min$(T).k$(K).out: $(FASTA) $(INTS)
	mkdir -p $(@D)
	$(SLIDER) \
		-l $(ELL) -d $(D) \
		-seq $(FASTA) -int $(INTS) \
		-o $@ \
		-m slider++ \
		-a $(K) \
		-st x2 \
		-min $(T)
	## remove dash in SLIDER's output
	tr '-' ' ' < $@.txt > $@
	rm -f $@.txt
