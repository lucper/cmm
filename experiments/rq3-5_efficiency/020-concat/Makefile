SHELL       := /usr/bin/bash
.SHELLFLAGS := -o pipefail -c

.PHONY: usage
usage:
	@echo "Usage: make all"
	@echo "       OUT=<output TSV>"
	@echo "       INPUTS=<.perf files from 010-measure>"

.DELETE_ON_ERROR:

.PHONY: all
all: $(OUT)

# Concatenates the measurements, prepending the trial number taken from each
# file name. The list of inputs can change without any input being newer than
# $(OUT) (e.g., when a node count is removed), so $(OUT) is always rebuilt, but
# it is only replaced when its contents change, so that the plots are not redone.
.PHONY: FORCE
$(OUT): $(INPUTS) FORCE
	awk '{trial = FILENAME; sub(/.*\.trial/, "", trial); sub(/\.perf$$/, "", trial); \
	      print trial "\t" $$0}' $(INPUTS) > $@.tmp
	if cmp -s $@.tmp $@; then rm $@.tmp; else mv $@.tmp $@; fi
