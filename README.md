# Code for "Mining Correlated Motifs with Wildcards in PPI networks: Exact and Optimal"

Here, you will find the source code for the exact algorithm proposed in "Mining Correlated Motifs with Wildcards in PPI networks: Exact and Optimal".
To reproduce Section 5 (**Experimental Results**), see [Reproducing the experiments](#reproducing-the-experiments).

## Requirements
- C++17 compiler GCC 7+
- GNU/Linux system (e.g., Ubuntu, Fedora)

## Installation
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

## Usage

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

## Example
The following command runs the program in the example input files above.
It retrieves the top 5 motif pairs with highest $\chi^2$-score.
Each motif will have length 5 and 2 wildcard symbols (`x`).
```bash
./cmm -s example.fa -i example.int -l 5 -d 2 -f x2 -k 5
```

This command yields the following output:
```text
MLHxx TQTxx 0,17
MLHxx QTFxx 0,17
MLHxx TFIxx 0,17
MLHxx TKDxx 0,17
MLHxx RTQxx 0,17
```

## Reproducing the experiments
The script `runme.sh` runs all the experiments of Section 5 (**Experimental Results**) and collects the results in `output/`.
It has two modes:
```bash
./runme.sh      # full experiments of the paper (takes several days with 128 threads)
./runme.sh -s   # subset evaluated by the ALENEX 2027 Artifact Evaluation Committee (takes a few hours with 64 threads)
```
All parameters of both modes are set at the top of `runme.sh`.
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
python3 python3-pip python3-venv
```

To run each experiment suite separately, see [`experiments/README.md`](experiments/README.md).

## Citation
If you use this software, please cite the following paper.
```text
Lucas Peres Oliveira, Ben Bals, Aalt-Jan van Dijk, and Solon P. Pissis. 
Mining Correlated Motifs with Wildcards in PPI networks: Exact and Optimal. 
In SIAM Symposium on Algorithm Engineering and Experiments, ALENEX 2027, 2027.
```
