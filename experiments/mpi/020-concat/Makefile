SHELL       := /usr/bin/bash
.SHELLFLAGS := -o pipefail -c

.PHONY: usage
usage:
	@echo "Usage: make all"
	@echo "       DATASET=<dataset name>"
	@echo "       L=<motif length>"
	@echo "       D=<number of wildcards>"
	@echo "       F=<support function>"
	@echo "       K=<number of top motif pairs>"

.DELETE_ON_ERROR:

NAME   := $(DATASET).layouts.l$(L)d$(D).$(F).k$(K)
OUT    := $(NAME).tsv
INPUTS := $(wildcard ../010-measure/$(DATASET)/$(NAME).*.perf)
PAIRS  := $(INPUTS:.perf=.out)

.PHONY: all
all: $(OUT) same-pairs

# Concatenates the measurements, prepending the nodes, the layout and the trial number taken from each
# file name. The list of inputs can change without any input being newer than $(OUT), so $(OUT) is always
# rebuilt, but it is only replaced when its contents change.
.PHONY: FORCE
$(OUT): $(INPUTS) FORCE
	@[ -n "$(INPUTS)" ] || { echo "FAILED: no measurements of $(NAME)"; exit 1; }
	awk '{n = split(FILENAME, parts, "."); nodes = parts[n - 3]; sub(/^n/, "", nodes); \
	      trial = parts[n - 1]; sub(/^trial/, "", trial); \
	      print nodes "\t" parts[n - 2] "\t" trial "\t" $$0}' $(INPUTS) > $@.tmp
	if cmp -s $@.tmp $@; then rm $@.tmp; else mv $@.tmp $@; fi

# Every run must find the same pairs, whatever its nodes and layout.
.PHONY: same-pairs
same-pairs:
	@for pairs in $(PAIRS); do \
	     cmp -s $(firstword $(PAIRS)) $$pairs || { echo "FAILED: $$pairs differs from $(firstword $(PAIRS))"; exit 1; }; \
	 done
