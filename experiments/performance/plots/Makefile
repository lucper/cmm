PLOT   := ../venv/bin/python3 plot.py

plots: rq1.pdf rq2.pdf rq3_time.pdf rq3_memory.pdf
.PHONY: plots

rq1.pdf: ../rq1.tsv
	$(PLOT) --input $< --rq rq1 --output $@

rq2.pdf: ../rq2.tsv
	$(PLOT) --input $< --rq rq2 --output $@

rq3_time.pdf rq3_memory.pdf: ../rq3.tsv
	$(PLOT) --input $< --rq rq3 --output rq3.pdf
