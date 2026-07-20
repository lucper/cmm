# Research questions

- **RQ1:** How do the scores of SLIDER's solution set compare to the ones of the optimal solution set?
- **RQ2:** What is the recall of SLIDER’s solution set relative to the optimal solution set?

## Prerequisites

Under the `experiments` directory:
1. Build the binaries by running `make -C src` (produces `bin/cmm_test`,
   `bin/cmm_eval`, `bin/cmm_dedup`).
2. Have the real datasets available under `data/real_instances/` (produced by
   `./data/run.sh`).

SLIDER is invoked from `010-soln-slider/` as a Java jar under
`000-slider/SliderLight/`. The plotting stages use `gnuplot` (with the
`cairolatex` terminal) and produce LaTeX/PGF `.tex` figures.

## Usage

```bash
./run.sh -D <dataset> -f <fasta> -i <ints> \
         -r <R> -t <T> -l <L> -d <D> \
         -k <K> -K <K_TOP> -H <h> -c <C> -p <P>
```

All options are required:

| Flag | Meaning |
|------|---------|
| `-D` | dataset name (used in every output filename) |
| `-f` | FASTA file with protein sequences |
| `-i` | interactions file (`u v` per line) |
| `-r` | number $R$ of runs for SLIDER |
| `-t` | max convergence time $T$ for SLIDER (minutes) |
| `-l` | motif length $\ell$ |
| `-d` | number of wildcards $d$ |
| `-k` | number $K$ of top motif pairs (per run) for SLIDER |
| `-K` | number of top motif pairs to keeo after aggregation |
| `-H` | similarity proximity $h$ |
| `-c` | coverage cutoff $C$ in $[0,100]$ |
| `-p` | threads for the exact algorithm |

The experiment is a sequence of numbered `make` stages, each consuming a previous stage's output.
`run.sh` drives them in order:

| Directory | Produces |
|-----------|----------|
| `010-soln-slider/` | SLIDER's solution sets |
| `010-soln-cmm/` | Exact algorithm's solution set |
| `020-agg/` | SLIDER solution set after aggregating top-$K$ motif pairs across $R$ runs |
| `025-dedup/` | SLIDER and exact algorithm's *reduced* solution sets |
| `030-score-curve/` | Normalized score curve for all methods |
| `030-eval/` | Table of highest similarity motif pairs between SLIDER and exact algorithm's reduced solution sets |
| `035-eval-density/` | Similarity density plot |
| `035-eval-coverage/` | Table with coverage of exact algorithm by SLIDER |

## Files

| File / dir | Role |
|------------|------|
| `run.sh` | Driver: runs all stages in order for one parameter cell. |
| `methods.csv` | Maps method keys (`m_slider`, `seq_slider`, `cmm`) to display names. |
| `<NNN>-*/` | The numbered pipeline stages, each with its own `Makefile`. |
| `000-slider/` | The SLIDER source code, invoked by `010-soln-slider`. |
| `<stage>/<dataset>/` | Per-stage outputs (generated). |

## Output format

Filenames encode the full parameter cell, so outputs from different cells never collide.
A representative name:

```
<dataset>.<method>.l<L>d<D>.agg.trials<R>.min<T>.k<K_TOP>.h<H>.dedup.tsv
```

where `<method>` is one of `seq_slider`, `m_slider`, or `cmm`.
The exact solution omits the `trials`/`min` fields, since it does not depend on $R$ or $T$:

```
<dataset>.cmm.l<L>d<D>.k<K_TOP>.out
```

## Plots

`030-score-curve` and `035-eval-density` each produce a figure plus a standalone legend, via `gnuplot` scripts (`plot_curve.gp`, `plot_density.gp`).
Method display names and colors are driven by `methods.csv`.
