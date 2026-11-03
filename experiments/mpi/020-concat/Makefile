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

# Concatenates the measurements, prepending the layout and the trial number taken from each file
# name. The list of inputs can change without any input being newer than $(OUT), so $(OUT) is always
# rebuilt, but it is only replaced when its contents change.
.PHONY: FORCE
$(OUT): $(INPUTS) FORCE
	awk '{n = split(FILENAME, parts, "."); trial = parts[n - 1]; sub(/^trial/, "", trial); \
	      print parts[n - 2] "\t" trial "\t" $$0}' $(INPUTS) > $@.tmp
	if cmp -s $@.tmp $@; then rm $@.tmp; else mv $@.tmp $@; fi
