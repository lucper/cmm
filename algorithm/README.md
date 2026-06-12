## Installation
**Requirements:**
- C++17 compiler GCC 7+
- GNU/Linux system (e.g., Ubuntu, Fedora)

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
