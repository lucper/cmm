PROG  = ../build/bpm
override INPUT = $(V) $(E) $(ell) $(d) $(k) $(mode)
ID ?= ""

all: profile$(ID).svg
.PHONY: all

clear:
	rm -rf *.data *.perf-folded *.svg *.out
.PHONY: clear

profile$(ID).svg: out$(ID).perf-folded
	./FlameGraph/flamegraph.pl $< > $@

out$(ID).perf-folded: perf$(ID).data
	perf script -i $< | ./FlameGraph/stackcollapse-perf.pl > $@

perf$(ID).data: FlameGraph
	perf record -F 99 -g -o $@ -- $(PROG) $(INPUT) > $(ID).out

FlameGraph:
	git clone https://github.com/brendangregg/$@.git
