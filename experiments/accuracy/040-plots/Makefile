SHELL := /usr/bin/bash

EVAL_DIR  := ../030-eval
AGG_DIR   := ../020-agg
DEDUP_DIR := ../025-dedup

## change this
EXACT_DIR := ../cmm_out

VENV_DIR := venv
PYTHON   := $(VENV_DIR)/bin/python
PIP      := $(VENV_DIR)/bin/pip

## default
## cells and datasets are defined on run.sh
DATASETS  := ecoli yeast fly plant human
T_VALUES  := 5 10 15 20 25 30
R_VALUES  := 5 10 15 20 25 30

## box plot cells
RQ1_R := 5
RQ1_T := 30
RQ2_R := 30
RQ2_T := 5

## eval_tsv(ds, R, T): non-dedup eval TSV.
eval_tsv       = $(EVAL_DIR)/$(1)/$(1).slider++.l$(ELL)d$(D).agg.trials$(2).min$(3).k$(K_TOP).h$(H).tsv
## eval_tsv_dedup(ds, R, T): dedup eval TSV.
eval_tsv_dedup = $(EVAL_DIR)/$(1)/$(1).slider++.l$(ELL)d$(D).agg.trials$(2).min$(3).k$(K_TOP).h$(H).dedup.tsv

## slider_agg(ds, R, T): raw slider aggregate (with X Y x2 header).
slider_agg   = $(AGG_DIR)/$(1)/$(1).slider++.l$(ELL)d$(D).agg.trials$(2).min$(3).k$(K_TOP).out
## slider_dedup(ds, R, T): deduped slider aggregate.
slider_dedup = $(DEDUP_DIR)/$(1)/$(1).slider++.l$(ELL)d$(D).agg.trials$(2).min$(3).k$(K_TOP).dedup.out

## exact_raw(ds): raw exact solution from cmm_out.
exact_raw   = $(EXACT_DIR)/$(1).cmm.l$(ELL)d$(D).k$(K_TOP).out
## exact_dedup(ds): deduped exact solution from 025-dedup.
exact_dedup = $(DEDUP_DIR)/$(1)/$(1).cmm.l$(ELL)d$(D).k$(K_TOP).dedup.out

RQ1_TSVS         = $(foreach d,$(DATASETS),$(foreach t,$(T_VALUES),$(call eval_tsv,$(d),$(RQ1_R),$(t))))
RQ1_TSVS_DEDUP   = $(foreach d,$(DATASETS),$(foreach t,$(T_VALUES),$(call eval_tsv_dedup,$(d),$(RQ1_R),$(t))))
RQ1_INPUTS       = $(foreach d,$(DATASETS),$(foreach t,$(T_VALUES),$(d):$(t):$(call eval_tsv,$(d),$(RQ1_R),$(t))))
RQ1_INPUTS_DEDUP = $(foreach d,$(DATASETS),$(foreach t,$(T_VALUES),$(d):$(t):$(call eval_tsv_dedup,$(d),$(RQ1_R),$(t))))

RQ2_TSVS         = $(foreach d,$(DATASETS),$(foreach r,$(R_VALUES),$(call eval_tsv,$(d),$(r),$(RQ2_T))))
RQ2_TSVS_DEDUP   = $(foreach d,$(DATASETS),$(foreach r,$(R_VALUES),$(call eval_tsv_dedup,$(d),$(r),$(RQ2_T))))
RQ2_INPUTS       = $(foreach d,$(DATASETS),$(foreach r,$(R_VALUES),$(d):$(r):$(call eval_tsv,$(d),$(r),$(RQ2_T))))
RQ2_INPUTS_DEDUP = $(foreach d,$(DATASETS),$(foreach r,$(R_VALUES),$(d):$(r):$(call eval_tsv_dedup,$(d),$(r),$(RQ2_T))))

RQ1_BOX_FILES         = $(foreach d,$(DATASETS),$(call slider_agg,$(d),$(RQ1_R),$(RQ1_T))   $(call exact_raw,$(d)))
RQ1_BOX_FILES_DEDUP   = $(foreach d,$(DATASETS),$(call slider_dedup,$(d),$(RQ1_R),$(RQ1_T)) $(call exact_dedup,$(d)))
RQ1_BOX_INPUTS        = $(foreach d,$(DATASETS),$(d):slider++:$(call slider_agg,$(d),$(RQ1_R),$(RQ1_T))   $(d):cmm:$(call exact_raw,$(d)))
RQ1_BOX_INPUTS_DEDUP  = $(foreach d,$(DATASETS),$(d):slider++:$(call slider_dedup,$(d),$(RQ1_R),$(RQ1_T)) $(d):cmm:$(call exact_dedup,$(d)))

RQ2_BOX_FILES         = $(foreach d,$(DATASETS),$(call slider_agg,$(d),$(RQ2_R),$(RQ2_T))   $(call exact_raw,$(d)))
RQ2_BOX_FILES_DEDUP   = $(foreach d,$(DATASETS),$(call slider_dedup,$(d),$(RQ2_R),$(RQ2_T)) $(call exact_dedup,$(d)))
RQ2_BOX_INPUTS        = $(foreach d,$(DATASETS),$(d):slider++:$(call slider_agg,$(d),$(RQ2_R),$(RQ2_T))   $(d):cmm:$(call exact_raw,$(d)))
RQ2_BOX_INPUTS_DEDUP  = $(foreach d,$(DATASETS),$(d):slider++:$(call slider_dedup,$(d),$(RQ2_R),$(RQ2_T)) $(d):cmm:$(call exact_dedup,$(d)))

.PHONY: usage all
usage:
	@echo "Usage: make DATASETS='ds1 ds2 ...' ELL=<l> D=<d> K_TOP=<k> H=<h>"
	@echo "            [THRESHOLD=<s>] <target>"
	@echo
	@echo "Targets:"
	@echo "  rq1_recall.h<H>.pdf         rq2_recall.h<H>.pdf"
	@echo "  rq1_recall.h<H>.dedup.pdf   rq2_recall.h<H>.dedup.pdf"
	@echo "  rq1_boxplot.pdf             rq2_boxplot.pdf"
	@echo "  rq1_boxplot.dedup.pdf       rq2_boxplot.dedup.pdf"
	@echo "  all"

.DELETE_ON_ERROR:

all: \
    rq1_recall.h$(H).thr$(THRESHOLD).pdf rq1_recall.h$(H).thr$(THRESHOLD).dedup.pdf \
    rq2_recall.h$(H).thr$(THRESHOLD).pdf rq2_recall.h$(H).thr$(THRESHOLD).dedup.pdf \
    rq1_boxplot.pdf      rq1_boxplot.dedup.pdf \
    rq2_boxplot.pdf      rq2_boxplot.dedup.pdf

rq1_recall.h$(H).thr$(THRESHOLD).pdf: $(RQ1_TSVS) | $(VENV_DIR)/bin/activate
	$(PYTHON) plot_top_n_recall.py \
		--input $(RQ1_INPUTS) \
		--output $@ \
		--threshold 0.$(THRESHOLD) \
		--x-label "Max convergence time (min)"

rq1_recall.h$(H).thr$(THRESHOLD).dedup.pdf: $(RQ1_TSVS_DEDUP) | $(VENV_DIR)/bin/activate
	$(PYTHON) plot_top_n_recall.py \
		--input $(RQ1_INPUTS_DEDUP) \
		--output $@ \
		--threshold 0.$(THRESHOLD) \
		--x-label "Max convergence time (min)"

rq2_recall.h$(H).thr$(THRESHOLD).pdf: $(RQ2_TSVS) | $(VENV_DIR)/bin/activate
	$(PYTHON) plot_top_n_recall.py \
		--input $(RQ2_INPUTS) \
		--output $@ \
		--threshold 0.$(THRESHOLD) \
		--x-label "Number of runs"

rq2_recall.h$(H).thr$(THRESHOLD).dedup.pdf: $(RQ2_TSVS_DEDUP) | $(VENV_DIR)/bin/activate
	$(PYTHON) plot_top_n_recall.py \
		--input $(RQ2_INPUTS_DEDUP) \
		--output $@ \
		--threshold 0.$(THRESHOLD) \
		--x-label "Number of runs"

rq1_boxplot.pdf: $(RQ1_BOX_FILES) | $(VENV_DIR)/bin/activate
	$(PYTHON) plot_score_dist_boxplot.py \
		--input $(RQ1_BOX_INPUTS) \
		--output $@ \
		--score-label 'Normalized $$f_{\chi^2}$$' \
		--baseline cmm

rq1_boxplot.dedup.pdf: $(RQ1_BOX_FILES_DEDUP) | $(VENV_DIR)/bin/activate
	$(PYTHON) plot_score_dist_boxplot.py \
		--input $(RQ1_BOX_INPUTS_DEDUP) \
		--output $@ \
		--score-label 'Normalized $$f_{\chi^2}$$' \
		--baseline cmm

rq2_boxplot.pdf: $(RQ2_BOX_FILES) | $(VENV_DIR)/bin/activate
	$(PYTHON) plot_score_dist_boxplot.py \
		--input $(RQ2_BOX_INPUTS) \
		--output $@ \
		--score-label 'Normalized $$f_{\chi^2}$$' \
		--baseline cmm

rq2_boxplot.dedup.pdf: $(RQ2_BOX_FILES_DEDUP) | $(VENV_DIR)/bin/activate
	$(PYTHON) plot_score_dist_boxplot.py \
		--input $(RQ2_BOX_INPUTS_DEDUP) \
		--output $@ \
		--score-label 'Normalized $$f_{\chi^2}$$' \
		--baseline cmm

$(VENV_DIR)/bin/activate: requirements.txt
	python3 -m venv $(VENV_DIR)
	$(PIP) install --upgrade pip
	$(PIP) install -r $<
