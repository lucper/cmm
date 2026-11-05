SHELL       := /usr/bin/bash
.SHELLFLAGS := -o pipefail -c

.PHONY: usage
usage:
	@echo "Usage: make setup, once, to install Julia and the packages of the script in .julia/ (Linux x86-64)"
	@echo "       make all"
	@echo "       DATASET=<dataset name>"
	@echo "       L=<motif length>"
	@echo "       D=<number of wildcards>"
	@echo "       F=<support function>"
	@echo "       K=<number of top motif pairs>"
	@echo "       TRIAL=<trial of the runs on 1 node to simulate from> (default: 1)"
	@echo "       JULIA=<another Julia program> (optional, instead of the one of make setup)"

.DELETE_ON_ERROR:

JULIA_VERSION := 1.13.0
JULIA_SHA256  := 8975da61c128a5e5ded3e719e868da8c8781deb7ad7913d37fb99be02a81904b
JULIA_URL     := https://julialang-s3.julialang.org/bin/linux/x64/1.13/julia-$(JULIA_VERSION)-linux-x86_64.tar.gz
JULIA         ?= .julia/julia-$(JULIA_VERSION)/bin/julia

# Everything Julia installs goes to .julia/, not to ~/.julia. The final colon keeps the libraries
# that come with Julia, already compiled. Project.toml and Manifest.toml fix the packages and versions.
export JULIA_DEPOT_PATH := $(CURDIR)/.julia/depot:
JULIA_RUN := $(JULIA) --project=$(CURDIR)

.PHONY: setup
setup: $(JULIA)
	$(JULIA_RUN) -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

$(JULIA):
	mkdir -p .julia
	curl -fL --retry 3 -o .julia/julia.tar.gz $(JULIA_URL)
	echo '$(JULIA_SHA256)  .julia/julia.tar.gz' | sha256sum -c -
	tar -xzf .julia/julia.tar.gz -C .julia
	rm .julia/julia.tar.gz
	touch $@

TRIAL ?= 1

NAME  := $(DATASET).layouts.l$(L)d$(D).$(F).k$(K)
OUT   := $(DATASET)/$(NAME).txt
TABLE := ../020-concat/$(NAME).tsv
# The times of each layout on 1 node. The script also reads the .perf next to each.
TIMES := $(wildcard ../010-measure/$(DATASET)/$(NAME).n1.*x*.trial$(TRIAL).tsv)

.PHONY: all
all: $(OUT)

# Writes the printed results to $(OUT), and for each layout <name>-<layout>.pdf and .png, and
# <name>-<layout>-speedup.pdf and .png, next to it. The trials are cut into one chunk per Julia
# thread, each with its own random numbers, so a fixed number of threads, with the same version of
# Julia, gives the same results.
$(OUT): simulate_tasks.jl $(TIMES) $(TIMES:.tsv=.perf) $(TABLE)
	@[ -n "$(TIMES)" ] || { echo "FAILED: no times of $(NAME) on 1 node, trial $(TRIAL)"; exit 1; }
	@[ -x "$(JULIA)" ] || { echo "FAILED: no Julia at $(JULIA): run make setup first"; exit 1; }
	mkdir -p $(@D)
	$(JULIA_RUN) --threads=6 simulate_tasks.jl $(TIMES) --results $(TABLE) --plot $(DATASET)/$(NAME) > $@.tmp
	mv $@.tmp $@
