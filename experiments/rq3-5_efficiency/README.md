## Research questions

- **RQ3:** How does our algorithm scale as $V$ (and thus $N$) grows?
- **RQ4:** How does our algorithm scale as $E$ grows?
- **RQ5:** How does our algorithm scale as the number $t$ of threads used increases?

## Prerequisites

1. Build `cmm_perf` by running `make`.
2. Run `../data/run.sh` to generate the artificial instances used for this experiment.

Python dependencies for plotting (`pandas`, `matplotlib`, `seaborn`) are listed
in `requirements.txt` and installed automatically into a local `venv/` on first
run.

## Usage

```bash
./run.sh              # run experiments, then plot
./run.sh --plot-only  # re-plot existing TSVs without re-running measurements
```

The run is resumable.
Each measurement is keyed by `(trial, instance, ell, d, threads)`; rows already present in the output TSV are skipped, so an interrupted run continues where it left off, and re-running only fills in what is missing.

## Files

| File                | Role                                                              |
|---------------------|-------------------------------------------------------------------|
| `run.sh`            | Driver: runs `cmm_perf` over the instances, then plots.           |
| `plot.py`           | Renders the TSVs into PDFs under `plots/`.                        |
| `requirements.txt`  | Python dependencies for `plot.py`.                                |
| `cmm_perf`          | The profiled binary version of `cmm` (exact algorithm).           |
| `rq{3,4,5}_out.tsv` | Measurement output, one row per run (generated).                  |
| `plots/`            | Rendered PDFs (generated).                                        |

## Output format

Each row of `rq{3,4,5}_out.tsv` is one `cmm_perf` run with the trial number
prepended, tab-separated:

| # | Column | | # | Column | | # | Column |
|---|--------|-|---|--------|-|---|--------|
| 1 | `trial` | | 6 | `ell` | | 11 | `pruning_cnt` |
| 2 | `instance` | | 7 | `d` | | 12 | `num_threads_requested` |
| 3 | `edge_density` | | 8 | `max_rank` | | 13 | `num_threads_spawned` |
| 4 | `V` | | 9 | `min_rank` | | 14 | `time_ms` |
| 5 | `N` | | 10 | `avg_rank` | | 15 | `peak_ram_kb` |

Note that `num_threads_spawned` (13) may be lower than `num_threads_requested` (12).

## Plots

`plot.py` writes PDFs into `plots/`, keyed by RQ. Time and memory are separate figures, with the legend saved as a standalone PDF.
The exact styling per RQ (log scaling, polynomial fits, `max_rank` annotations, RQ5 speedup labels) is set by the `plot.py` invocations in `run.sh`.
