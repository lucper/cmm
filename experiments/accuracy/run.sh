#!/usr/bin/bash

set -euo pipefail

DATASET=$1
FASTA=$(realpath $2)
INTS=$(realpath $3)

L=8
D=3
K=10000
K_TOP=10000
# set to 0 for stricter evaluation
H=8
C=70

## (trials, time) cells: RQ1 (vary time, trials=5) + RQ2 (vary trials, time=5).
EXP_CELLS=(
    #"5:5"  "5:10" "5:15" "5:20" "5:25" "5:30"
    #"10:5" "15:5" "20:5" "25:5" "30:5"
    "10:10"
)

## Compute max trials needed at each time value (so RQ2's t=5 column doesn't
## redundantly trigger SLIDER runs at every other time).
declare -A MAX_R_AT_T
for cell in "${EXP_CELLS[@]}"
do
    r=${cell%:*}; t=${cell#*:}
    if [ -z "${MAX_R_AT_T[$t]:-}" ] || [ "${MAX_R_AT_T[$t]}" -lt "${r}" ]
    then
        MAX_R_AT_T[$t]=$r
    fi
done

echo "=== Stage 010-soln: SLIDER runs ===" >&2
for t in "${!MAX_R_AT_T[@]}"
do
    for i in $(seq 1 ${MAX_R_AT_T[$t]})
    do
        make -C 010-soln all DATASET=${DATASET} FASTA=${FASTA} INTS=${INTS} I=${i} T=${t} L=${L} D=${D} K=${K}
    done
done

echo "=== Stage 020-agg: Aggregation ===" >&2
for cell in "${EXP_CELLS[@]}"
do
    r=${cell%:*}; t=${cell#*:}
    make -C 020-agg all DATASET=${DATASET} R=${r} T=${t} L=${L} D=${D} K=${K} K_TOP=${K_TOP} -n
done

echo "=== Stage 025-dedup: Deduplication ===" >&2
## dedup SLIDER
for cell in "${EXP_CELLS[@]}"
do
    r=${cell%:*}; t=${cell#*:}
    make -C 025-dedup slider DATASET=${DATASET} FASTA=${FASTA} INTS=${INTS} R=${r} T=${t} L=${L} D=${D} K_TOP=${K_TOP} -n
done
## dedup Exact
make -C 025-dedup cmm DATASET=${DATASET} L=${L} D=${D} K_TOP=${K_TOP} -n


echo "=== Stage 030-score: Score evaluation ===" >&2
for cell in "${EXP_CELLS[@]}"
do
    r=${cell%:*}; t=${cell#*:}

    make -C 030-score all DATASET=${DATASET} R=${r} T=${t} L=${L} D=${D} K_TOP=${K_TOP} -n
done

echo "=== Stage 030-eval: Similarity evaluation ===" >&2
for cell in "${EXP_CELLS[@]}"
do
    r=${cell%:*}; t=${cell#*:}

    make -C 030-eval all DATASET=${DATASET} FASTA=${FASTA} INTS=${INTS} R=${r} T=${t} L=${L} D=${D} K_TOP=${K_TOP} H=${H} -n
done

echo "=== Stage 035-eval-recall: Recall against exact ===" >&2
for cell in "${EXP_CELLS[@]}"
do
    r=${cell%:*}; t=${cell#*:}

    make -C 035-eval-recall all DATASET=${DATASET} R=${r} T=${t} L=${L} D=${D} K_TOP=${K_TOP} H=${H} C=${C} -n
done

#echo "=== Stage 040-plots: Plotting ===" >&2
#make -C 040-plots \
#    DATASETS="${DATASETS[*]}" \
#    L=${L} D=${D} K_TOP=${K_TOP} H=${H} THRESHOLD=${THRESHOLD} \
#    all
