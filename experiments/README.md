# Experiments

Here, you will find two experiment suites to answer the five research questions stated in Section 5 (**Experimental Results**) of the paper "Mining Correlated Motifs with Wildcards in PPI networks: Exact and Optimal":

- Suite 1: `rq1-2_accuracy`
    - **RQ1:** How do the scores of SLIDER's solution set compare to the ones of the optimal solution set?
    - **RQ2:** What is the recall of SLIDER’s solution set relative to the optimal solution set?
- Suite 2: `rq3-5_efficiency`
    - **RQ3:** How does our algorithm scale as $V$ (and thus $N$) grows?
    - **RQ4:** How does our algorithm scale as $E$ grows?
    - **RQ5:** How does our algorithm scale as the number $t$ of threads used increases?

## One-time setup

Developed and run on a GNU/Linux system (x86-64).
The scripts use GNU tooling and bash features, so a Linux environment (or WSL) is assumed.

From this directory:

```bash
make -C src            # builds bin/cmm_test, bin/cmm_eval, bin/cmm_dedup
./data/run.sh          # fetches + cleans real data and generates artificial instances
```

Both suites call binaries in `bin/` and expect their inputs under `data/`.
The script `data/run.sh` populates both `data/real_instances/` (used by RQ1-RQ2) and `data/artificial_instances/` (used by RQ3-RQ5).

## RQ1 and RQ2 (accuracy experiment)

Each `run.sh` invocation runs the experiment for **one parameter cell** on a real dataset.
For example, from `rq1-2_accuracy`:
1. Pick a cleaned dataset from `data/real_instances/` (here taxon `9606`).
```bash
FA=../data/real_instances/9606.sequences.s700.clean.fa
INT=../data/real_instances/9606.links.s700.clean.int
```
2. Run a parameter cell:
```bash
./run.sh -D 9606 -f "$FA" -i "$INT" \
         -r 5 -t 10 -l 8 -d 3 -k 1000 -K 1000 -H 8 -c 70 -p 8
```

The `-l/-d` (motif length / number of wildcards) and `-k/-K/-H/-c` values above are illustrative
See `rq1-2_accuracy/README.md` for what every flag means and how the outputs are named.

## RQ3, RQ4, and RQ5 (efficiency experiment)

A single driver runs the exact algorithm over the artificial instances and plots the results:

```bash
cd rq3-5_efficiency
./run.sh               # run code for RQ3, RQ4, RQ5, then plot
./run.sh --plot-only   # re-plot existing results without re-measuring
```

The run is resumable: completed measurements are skipped, so an interrupted run continues where it left off.
Outputs land in `rq3-5_efficiency/plots/`.
See `rq3-5_efficiency/README.md` for the parameters and output format.
