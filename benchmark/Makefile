PROG  = ../build/bpm
INPUT = $(V) $(E) $(ell) $(d) $(k)

all: profile.svg
.PHONY: all

clear:
	rm -rf *.data *.perf-folded *.svg
.PHONY: clear

profile.svg: out.perf-folded
	./FlameGraph/flamegraph.pl $< > $@

out.perf-folded: perf.data
	perf script | ./FlameGraph/stackcollapse-perf.pl > $@

perf.data: FlameGraph
	perf record -F 99 -g -- $(PROG) $(INPUT)

FlameGraph:
	git clone https://github.com/brendangregg/$@.git
