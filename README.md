# Code for "Mining Correlated Motifs with Wildcards in PPI Networks: Exact and Optimal"

Here, you will find the source code for the exact algorithm proposed in "Mining Correlated Motifs with Wildcards in PPI Networks: Exact and Optimal".
To reproduce Section 5 (**Experimental Results**), see [Reproducing the experiments](#reproducing-the-experiments).

## Users' Guide

### Requirements
- C++17 compiler GCC 7+
- GNU/Linux system (e.g., Ubuntu, Fedora)
- For the [MPI version](#mpi-version) only: an MPI implementation with `mpicxx` (e.g., Open MPI)

### Installation
Run the following command in the current directory to compile and link the code:
```bash
make
```

The result will be the executable `cmm`.
By running `./cmm -h`, the followling message should be displayed.
```text
Correlated Motif Miner

Usage:
cmm [OPTION...]

  -s, --sequences arg           Path to FASTA file with protein sequences.
                                [required]
  -i, --interactions arg        Path to text file with lines formatted as 'u v',
                                where u and v are sequence IDs from the FASTA
                                file. [required]
  -l, --motif-length arg        Motif length. [required]
  -d, --number-of-wildcards arg
                                Number in [0,l) of wildcards in motif.
                                [required]
  -f, --support-function arg    Support function to sort motifs ('E', 'x2').
                                [required]
  -k, --number-of-motifs arg    Number of top k motifs. (default: 1)
  -t, --threads arg             Number of threads. (default: 1)
  -v, --version                 Print version.
  -h, --help                    Print usage.
```

### Usage

- Option `-s`:
FASTA file containing protein sequences.
For example:
```text
>seqA
MLHCARRYMLVRPRL
EYL
>seqB
MDRTQTFIKDCLFTK
DSALRYSILG
>seqC
MELSSPSKKTTTSPI
SSKS
```

- Option `-i`:
Text file containing protein interactions.
Specifically, it is expected a two-column space-separated text file in which each line contains a pair of sequence identifiers from the FASTA file provided to option `-s`.
For example:
```text
seqA seqB
seqC seqA
```

- Options `-l` and `d` are, respectively, the motif length and the number of wildcard symbols in the motif.

- Option `-f`:
Support function to be maximized by the algorithm.
We support `E`, the number of edges in which a motif pair co-occurs, and `x2`, the $\chi^2$-score.

### Example
The following command runs the program in the example input files above.
It retrieves the top 5 motif pairs with highest $\chi^2$-score.
Each motif will have length 5 and 2 wildcard symbols (`x`).
```bash
./cmm -s example.fa -i example.int -l 5 -d 2 -f x2 -k 5
```

This command yields the following output:
```text
ALRxx ARRxx 0.17
ALRxx ARxYx 0.17
ALRxx ARxxM 0.17
ALRxx AxRYx 0.17
ALRxx AxRxM 0.17
```

Motif pairs with equal scores are listed in lexicographic order of the first motif, then of the second (byte order, so the wildcard `x` comes after the letters).

### MPI version
`cmm_mpi` runs the search on several processes, possibly on several nodes, each with its own threads.
It is built separately, so `cmm` does not need MPI:
```bash
make cmm_mpi
```

It takes the same options as `cmm` and prints the same motif pairs, for any number of processes.
Option `-t` is the number of threads **per process**.
Only the first process reads the input files, so the other nodes do not need them.
For example, on 2 processes with 4 threads each:
```bash
mpirun -np 2 --bind-to none ./cmm_mpi -s example.fa -i example.int -l 5 -d 2 -f x2 -k 5 -t 4
```
`--bind-to none` matters on one machine: with one or two processes, Open MPI may otherwise bind each process to a single core, which all its threads would share.
On a cluster, bind each process to its own cores instead, as the job scheduler allows; for Slurm, see [`experiments/mpi/check_snellius.sh`](experiments/mpi/check_snellius.sh).

Each process searches its own share of the wildcard configurations, so there is little to gain from more processes than configurations: a motif length of 5 has at most 55, and a length of 8 has up to 2485.

The tests of the MPI version run with `make test_mpi`.

## Reproducing the experiments
The script `runme.sh` runs all the experiments of Section 5 (**Experimental Results**) and collects the results in `output/full/` (or `output/subset/` with `-s`).
It has two modes:
```bash
./runme.sh      # full experiments of the paper
./runme.sh -s   # subset for evaluation purposes
```
All parameters of both modes are set at the top of `runme.sh`.
The results are collected in `output/full/` (or `output/subset/` with `-s`), where `<L>` is the motif length and `<D>` the number of wildcards:

| File | Content | Paper |
|------|---------|-------|
| `rq1_score_curves_l<L>d<D>.pdf` | Normalized score curves of SLIDER and the exact algorithm, one panel per dataset (RQ1) | Figure 5.1 |
| `rq2_similarity_densities_l<L>d<D>.pdf` | Densities of the similarity between SLIDER's and the exact solution sets, one panel per dataset (RQ2) | Figure 5.2 |
| `rq2_coverage.tsv` | Coverage of the exact solution set by SLIDER, per dataset, $D$, and SLIDER variant (RQ2) | Table C.2 |
| `exact_algorithm_performance.tsv` | Running time and peak memory (total and per thread) of the exact algorithm, per dataset and $D$ | Table C.3 |
| `rq3-5_efficiency.pdf` | Running time and peak memory of the exact algorithm on the artificial instances (RQ3, RQ4, RQ5) | Figures 5.3, 5.4, and 5.5 |

The subset (`-s`) produces a selection of the rows and panels of these figures and tables, with smaller parameters for RQ3-RQ5.

Internet access is needed on the first run, since the datasets are downloaded from the STRING database.
On a Linux host, the run can be resumed: rerunning `runme.sh` skips the work that is already done (in Docker with `--rm`, intermediate results are discarded when the container exits).

To run `runme.sh` in Docker, use the provided `Dockerfile`:
```bash
docker build -t cmm-repro .
docker run --rm -it -v "$PWD/output:/work/output" cmm-repro               # full experiments
docker run --rm -it -v "$PWD/output:/work/output" cmm-repro ./runme.sh -s # subset
```

To run it directly on a Linux host, install the following (Ubuntu package names):
```text
build-essential make wget curl gzip gawk tar unzip
default-jre-headless gnuplot-nox
texlive-latex-base texlive-pictures texlive-latex-recommended
python3
```

To run each experiment suite separately, see [`experiments/README.md`](experiments/README.md).

## Citation
If you use this software, please cite the following paper.
```bibtex
@inproceedings{OBvDP27,
  author    = {Lucas P. Oliveira and Ben Bals and Aalt-Jan van Dijk and Solon P. Pissis},
  title     = {Mining Correlated Motifs with Wildcards in {PPI} Networks: Exact and Optimal},
  booktitle = {SIAM Symposium on Algorithm Engineering and Experiments, {ALENEX} 2027},
  publisher = {{SIAM}},
  year      = {2027}
}
```

The software is archived on Zenodo:
```bibtex
@software{cmm,
  author    = {Lucas P. Oliveira and Ben Bals and Aalt-Jan van Dijk and Solon P. Pissis},
  title     = {Mining Correlated Motifs with Wildcards in {PPI} Networks: Exact and Optimal},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.23017579},
  url       = {https://doi.org/10.5281/zenodo.23017579}
}
```
