SHELL   := /usr/bin/bash

SLIDER  := java -cp ../000-slider/SliderLight/dist/SliderLight.jar Framework.Framework

## change this
DATA    := ../../string_data_cleaned
FASTA   := $(DATA)/$(DATASET).s700.cleaned.fa
INTS    := $(DATA)/$(DATASET).s700.cleaned.int

.PHONY: usage
usage:
	@echo "Usage: make all"
	@echo "       DATASET=<dataset name>"
	@echo "       I=<run index>"
	@echo "       T=<max convergence time>"
	@echo "       L=<motif length>"
	@echo "       D=<number of wildcards>"
	@echo "       K=<number of top motif pairs>"

.DELETE_ON_ERROR:

.PHONY: all
all: $(DATASET)/$(DATASET).seq_slider.l$(L)d$(D).trial$(I).min$(T).k$(K).out \
	 $(DATASET)/$(DATASET).m_slider.l$(L)d$(D).trial$(I).min$(T).k$(K).out

$(DATASET)/$(DATASET).seq_slider.l$(L)d$(D).trial$(I).min$(T).k$(K).out: $(FASTA) $(INTS)
	mkdir -p $(@D)
	$(SLIDER) \
		-l $(L) -d $(D) \
		-seq $(FASTA) -int $(INTS) \
		-o $@ -a $(K) -st x2 -min $(T) \
		-m slider++
	## remove dash in SLIDER's output
	tr '-' ' ' < $@.txt > $@
	rm -f $@.txt

$(DATASET)/$(DATASET).m_slider.l$(L)d$(D).trial$(I).min$(T).k$(K).out: $(FASTA) $(INTS)
	mkdir -p $(@D)
	$(SLIDER) \
		-l $(L) -d $(D) \
		-seq $(FASTA) -int $(INTS) \
		-o $@ -a $(K) -st x2 -min $(T) \
		-m HCTimed \
		-f Slide
	## remove dash in SLIDER's output
	tr '-' ' ' < $@.txt > $@
	rm -f $@.txt
