SHELL   := /usr/bin/bash

SLIDER  := java -cp ../000-slider/SliderLight/dist/SliderLight.jar Framework.Framework

## change this
DATA    := ../../string_data_cleaned
FASTA   := $(DATA)/$(DATASET).s700.cleaned.fa
INTS    := $(DATA)/$(DATASET).s700.cleaned.int

ELL     := 8
D       := 3
K       := 10000

.DELETE_ON_ERROR:

$(DATASET)/$(DATASET).slider++.trial$(TRIAL).min$(T).k$(K).out: $(FASTA) $(INTS)
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
