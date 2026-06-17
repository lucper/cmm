#!/usr/bin/bash

set -euo pipefail

DATASET=$1
FASTA=$(realpath $2)
INTS=$(realpath $3)
L=$4
D=$5
K=$6
K_TOP=$7
R=$8 # set to 0 for stricter evaluation
T=$9
H=$10
C=$11

echo "=== Stage 010-soln: SLIDER runs ===" >&2
for i in $(seq 1 ${R})
do
    make -C 010-soln all \
        DATASET=${DATASET} FASTA=${FASTA} INTS=${INTS} I=${i} T=${T} L=${L} D=${D} K=${K}
done

echo "=== Stage 020-agg: Aggregation ===" >&2
make -C 020-agg all \
    DATASET=${DATASET} R=${R} T=${T} L=${L} D=${D} K=${K} K_TOP=${K_TOP}

echo "=== Stage 025-dedup: Deduplication ===" >&2
make -C 025-dedup slider \
    DATASET=${DATASET} FASTA=${FASTA} INTS=${INTS} R=${R} T=${T} L=${L} D=${D} K_TOP=${K_TOP}
make -C 025-dedup cmm \
    DATASET=${DATASET} FASTA=${FASTA} INTS=${INTS} L=${L} D=${D} K_TOP=${K_TOP}

echo "=== Stage 030-score: Score evaluation ===" >&2
make -C 030-score all \
    DATASET=${DATASET} R=${R} T=${T} L=${L} D=${D} K_TOP=${K_TOP}

echo "=== Stage 030-eval: Similarity evaluation ===" >&2
make -C 030-eval all \
    DATASET=${DATASET} FASTA=${FASTA} INTS=${INTS} R=${R} T=${T} L=${L} D=${D} K_TOP=${K_TOP} H=${H}

echo "=== Stage 035-eval-recall: Recall against exact ===" >&2
make -C 035-eval-recall all \
    DATASET=${DATASET} R=${R} T=${T} L=${L} D=${D} K_TOP=${K_TOP} H=${H} C=${C}

#echo "=== Stage 040-plots: Plotting ===" >&2
#make -C 040-plots \
#    DATASETS="${DATASETS[*]}" \
#    L=${L} D=${D} K_TOP=${K_TOP} H=${H} THRESHOLD=${THRESHOLD} \
#    all
