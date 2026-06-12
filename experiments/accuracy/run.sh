#!/usr/bin/bash

set -euo pipefail

DATASETS=(ecoli yeast fly plant human)

ELL=8
D=3
K_SLIDER=10000
K_TOP=1000
# set to 0 for stricter evaluation
H=0
THRESHOLD=70

## (trials, time) cells: RQ1 (vary time, trials=5) + RQ2 (vary trials, time=5).
EXP_CELLS=(
    "5:5"  "5:10" "5:15" "5:20" "5:25" "5:30"
    "10:5" "15:5" "20:5" "25:5" "30:5"
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
for ds in "${DATASETS[@]}"
do
    for t in "${!MAX_R_AT_T[@]}"
    do
        for i in $(seq 1 ${MAX_R_AT_T[$t]})
        do
            make -C 010-soln \
                DATASET=${ds} I=${i} T=${t} ELL=${ELL} D=${D} K=${K_SLIDER} \
                "${ds}/${ds}.slider++.l${ELL}d${D}.trial${i}.min${t}.k${K_SLIDER}.out"
        done
    done
done

echo "=== Stage 020-agg: Aggregation ===" >&2
for ds in "${DATASETS[@]}"
do
    for cell in "${EXP_CELLS[@]}"
    do
        r=${cell%:*}; t=${cell#*:}
        make -C 020-agg \
            DATASET=${ds} R=${r} T=${t} ELL=${ELL} D=${D} K=${K_SLIDER} K_TOP=${K_TOP} \
            "${ds}/${ds}.slider++.l${ELL}d${D}.agg.trials${r}.min${t}.k${K_TOP}.out"
    done
done

echo "=== Stage 025-dedup: Deduplication ===" >&2
## dedup SLIDER
for ds in "${DATASETS[@]}"
do
    for cell in "${EXP_CELLS[@]}"
    do
        r=${cell%:*}; t=${cell#*:}
        make -C 025-dedup \
            DATASET=${ds} R=${r} T=${t} ELL=${ELL} D=${D} K_TOP=${K_TOP} \
            "${ds}/${ds}.slider++.l${ELL}d${D}.agg.trials${r}.min${t}.k${K_TOP}.dedup.out"
    done
done
## dedup Exact
for ds in "${DATASETS[@]}"
do
    make -C 025-dedup \
        DATASET=${ds} ELL=${ELL} D=${D} K_TOP=${K_TOP} \
        "${ds}/${ds}.cmm.l${ELL}d${D}.k${K_TOP}.dedup.out"
done

echo "=== Stage 030-eval: Similarity evaluation ===" >&2
for ds in "${DATASETS[@]}"
do
    for cell in "${EXP_CELLS[@]}"
    do
        r=${cell%:*}; t=${cell#*:}

        # no dedup
        exact_sol="../cmm_out/${ds}.cmm.l${ELL}d${D}.k${K_TOP}.out"
        make -C 030-eval \
            DATASET=${ds} R=${r} T=${t} ELL=${ELL} D=${D} K_TOP=${K_TOP} H=${H} EXACT_SOL=${exact_sol} \
            "${ds}/${ds}.slider++.l${ELL}d${D}.agg.trials${r}.min${t}.k${K_TOP}.h${H}.tsv"

        # dedup
        make -C 030-eval \
            DATASET=${ds} R=${r} T=${t} ELL=${ELL} D=${D} K_TOP=${K_TOP} H=${H} \
            "${ds}/${ds}.slider++.l${ELL}d${D}.agg.trials${r}.min${t}.k${K_TOP}.h${H}.dedup.tsv"
    done
done

echo "=== Stage 035-recall: Recall against exact ===" >&2
## 035-recall consumes 030-eval at h=0 regardless of the H used in stage
## 030-eval above. Force the h=0 dedup tsv to exist before recall.
for ds in "${DATASETS[@]}"
do
    for cell in "${EXP_CELLS[@]}"
    do
        r=${cell%:*}; t=${cell#*:}

        make -C 030-eval \
            DATASET=${ds} R=${r} T=${t} ELL=${ELL} D=${D} K_TOP=${K_TOP} H=0 \
            "${ds}/${ds}.slider++.l${ELL}d${D}.agg.trials${r}.min${t}.k${K_TOP}.h0.dedup.tsv"

        make -C 035-recall \
            DATASET=${ds} R=${r} T=${t} ELL=${ELL} D=${D} K_TOP=${K_TOP} \
            "${ds}/${ds}.recall.l${ELL}d${D}.trials${r}.min${t}.tsv"
    done
done

echo "=== Stage 040-plots: Plotting ===" >&2
make -C 040-plots \
    DATASETS="${DATASETS[*]}" \
    ELL=${ELL} D=${D} K_TOP=${K_TOP} H=${H} THRESHOLD=${THRESHOLD} \
    all
