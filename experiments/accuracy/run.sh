#!/usr/bin/bash

set -euo pipefail

DATASETS=(ecoli yeast fly plant) # add human
ELL=8
D=3
K_SLIDER=10000
K_TOP=1000
H_VALUES=(0 ${ELL})

DEDUP=${DEDUP:-0}                # off by default, override at the call site with: DEDUP=1 ./run.sh
CMM_DIR="cmm_out"                # change this
OUT_DIR="slider_out"             # change this

CMM_EVAL=./cmm_eval
CMM_DEDUP=./cmm_dedup

agg_and_eval() {
    local ds=$1
    local t=$2
    local r=$3

    local fa="${SEQ_DIR}/${ds}.s700.cleaned.fa"
    local int="${SEQ_DIR}/${ds}.s700.cleaned.int"
    local exact="${CMM_DIR}/${ds}.cmm.out"

    ## Merge+sort top K of r trials at time t.
    local agg="${OUT_DIR}/${ds}/agg/${ds}.slider++.agg.min${t}.trials${r}.k${K_TOP}.out"
    if [ ! -f "${agg}" ]
    then
        echo "Aggregating ${ds} min${t} trials${r}" >&2
        local trial_files=()
        for i in $(seq 1 ${r})
        do
            trial_files+=("${OUT_DIR}/${ds}/${ds}.slider++.trial${i}.min${t}.k${K_SLIDER}.out")
        done
        cat "${trial_files[@]}" | sort -k3 -gr | awk -v n=${K_TOP} 'NR<=n' > "${agg}.tmp"
        mv "${agg}.tmp" "${agg}"
    else
        echo "Skipping agg ${ds} min${t} trials${r} -- already done" >&2
    fi

    ## No-dedup eval at each h.
    for h in "${H_VALUES[@]}"
    do
        local tsv="${OUT_DIR}/${ds}/eval/${ds}.slider++.agg.min${t}.trials${r}.k${K_TOP}.h${h}.tsv"
        if [ ! -f "${tsv}" ]
        then
            echo "Evaluating ${ds} min${t} trials${r} h${h} no-dedup" >&2
            ${CMM_EVAL} "${fa}" "${int}" \
                "${agg}" "${exact}" \
                ${h} "${tsv}.tmp" /dev/null
            mv "${tsv}.tmp" "${tsv}"
        else
            echo "Skipping eval ${ds} min${t} trials${r} h${h} no-dedup -- already done" >&2
        fi
    done

    ## Optional: dedup the top K and eval the dedup'd version at each h.
    if [ ${DEDUP} -ne 1 ]
    then
        return
    fi

    local agg_dedup="${OUT_DIR}/${ds}/agg/${ds}.slider++.agg.min${t}.trials${r}.k${K_TOP}.dedup.out"
    if [ ! -f "${agg_dedup}" ]
    then
        echo "Deduping ${ds} min${t} trials${r}" >&2
        ${CMM_DEDUP} "${fa}" "${int}" "${agg}" 0 "${agg_dedup}.tmp"
        mv "${agg_dedup}.tmp" "${agg_dedup}"
    else
        echo "Skipping dedup ${ds} min${t} trials${r} -- already done" >&2
    fi

    for h in "${H_VALUES[@]}"
    do
        local tsv="${OUT_DIR}/${ds}/eval/${ds}.slider++.agg.min${t}.trials${r}.k${K_TOP}.h${h}.dedup.tsv"
        if [ ! -f "${tsv}" ]
        then
            echo "Evaluating ${ds} min${t} trials${r} h${h} dedup" >&2
            ${CMM_EVAL} "${fa}" "${int}" \
                "${agg_dedup}" "${exact}" \
                ${h} "${tsv}.tmp" /dev/null
            mv "${tsv}.tmp" "${tsv}"
        else
            echo "Skipping eval ${ds} min${t} trials${r} h${h} dedup -- already done" >&2
        fi
    done
}

## RQ1: vary max time (trials fixed at 5)
echo "Running experiment for RQ1..." >&2
TRIALS=5
TIMES=(5 10 15 20 25 30)
for ds in "${DATASETS[@]}"
do
    for trial in $(seq 1 ${TRIALS})
    do
        for t in "${TIMES[@]}"
        do
            make -C 010-soln DATASET=${ds} TRIAL=${trial} T=${t} ELL=${ELL} D=${D}
        done
    done
    for t in "${TIMES[@]}"
    do
        ## agg_and_eval ${ds} ${t} ${TRIALS}
    done
done

## RQ2: vary number of trials (max time fixed at 5)
#echo "Running experiment for RQ2..." >&2
#TIME=5
#TRIAL_COUNTS=(5 10 15 20 25 30)
#MAX_TRIALS=30
#for ds in "${DATASETS[@]}"
#do
#    for trial in $(seq 1 ${MAX_TRIALS})
#    do
#        run_slider ${ds} ${trial} ${TIME}
#    done
#    for r in "${TRIAL_COUNTS[@]}"
#    do
#        agg_and_eval ${ds} ${TIME} ${r}
#    done
#done

echo "All experiments done." >&2
