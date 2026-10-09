SHELL       := /usr/bin/bash
.SHELLFLAGS := -o pipefail -c

INPUT      := ../020-concat
PLOT       := ./plot.py
ANNOTATE   := 8:4 8:7 5:2 5:4   # points annotated with max_rank
POLY_PAIRS := 5:1 5:3 8:2 8:5   # points given a polynomial fit on RQ3

.PHONY: usage
usage:
	@echo "Usage: make all"
	@echo "       TSV=<TSV name in 020-concat (rq3.*, rq4.*, or rq5.*)>"
	@echo "       RQ5_ELL=<motif length for RQ5>"

.DELETE_ON_ERROR:

# plot.py writes several PDFs per TSV, named <TSV name without .tsv>_*.pdf; the
# .plotted file marks that they are up to date. --output is given with a .pdf
# extension because plot.py strips the last extension and the names contain dots.
.PHONY: all
all: $(TSV:.tsv=.plotted)

# RQ3: log-scale both axes, quadratic fit on time, linear fit on memory per thread.
rq3.%.plotted: $(INPUT)/rq3.%.tsv $(PLOT)
	$(PLOT) \
		--input $< \
		--rq rq3 \
		--output rq3.$*.pdf \
		--log-scale-time \
		--log-scale-mem \
		--fit-poly-time $(POLY_PAIRS) \
		--fit-poly-mem $(POLY_PAIRS) \
		--poly-degree-time 2 \
		--poly-degree-mem 1 \
		--no-poly-label \
		--max-rank $(ANNOTATE) \
		--mem-per-thread
	touch $@

# RQ4: annotate selected points with max_rank (linear axes, no fits), memory per thread.
rq4.%.plotted: $(INPUT)/rq4.%.tsv $(PLOT)
	$(PLOT) \
		--input $< \
		--rq rq4 \
		--output rq4.$*.pdf \
		--max-rank $(ANNOTATE) \
		--mem-per-thread
	touch $@

# RQ5: speedup labels, log-scale both axes, total memory.
rq5.%.plotted: $(INPUT)/rq5.%.tsv $(PLOT)
	$(PLOT) \
		--input $< \
		--rq rq5 \
		--output rq5.$*.pdf \
		--speedup \
		--rq5-ell $(RQ5_ELL) \
		--log-scale-time \
		--log-scale-mem
	touch $@
