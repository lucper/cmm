#!/usr/bin/bash

set -euo pipefail

DATASETS=(ecoli yeast fly plant) # add human
ELL=8
D=3
K_SLIDER=10000
K_TOP=1000
H_VALUES=(0 ${ELL})

SEQ_DIR="../string_data_cleaned" # change this
CMM_DIR="cmm_out"                # change this
OUT_DIR="slider_out"             # change this

CMM_EVAL=./cmm_eval
SLIDER_JAR="../competitors/SliderLight/dist/SliderLight.jar"

## Sanity checks
[ -x "${CMM_EVAL}" ] || { echo "Error: ${CMM_EVAL} not found or not executable" >&2; exit 1; }
[ -f "${SLIDER_JAR}" ] || { echo "Error: SliderLight jar not found at ${SLIDER_JAR}" >&2; exit 1; }

for ds in "${DATASETS[@]}"
do
    for ext in fa int
    do
        f="${SEQ_DIR}/${ds}.s700.cleaned.${ext}"
        [ -f "${f}" ] || { echo "Error: ${f} not found" >&2; exit 1; }
    done
    f="${CMM_DIR}/${ds}.cmm.out"
    [ -f "${f}" ] || { echo "Error: ${f} (EXACT solution) not found" >&2; exit 1; }
done

for ds in "${DATASETS[@]}"
do
    mkdir -p "${OUT_DIR}/${ds}" "${OUT_DIR}/${ds}/agg" "${OUT_DIR}/${ds}/eval"
done

run_slider() {
    local ds=$1
    local trial=$2
    local t=$3

    local out="${OUT_DIR}/${ds}/${ds}.slider++.trial${trial}.min${t}.k${K_SLIDER}.out"
    if [ -f "${out}" ]
    then
        echo "Skipping SLIDER ${ds} trial${trial} min${t} -- already done" >&2
        return
    fi

    echo "Running SLIDER ${ds} trial${trial} min${t}" >&2
    local fa="${SEQ_DIR}/${ds}.s700.cleaned.fa"
    local int="${SEQ_DIR}/${ds}.s700.cleaned.int"
    local raw="${out}.raw"
    java -cp "${SLIDER_JAR}" Framework.Framework \
        -l ${ELL} -d ${D} \
        -seq "${fa}" \
        -int "${int}" \
        -o "${raw}" \
        -m slider++ -a ${K_SLIDER} -st x2 -min ${t}
    tr '-' ' ' < "${raw}.txt" > "${out}.tmp"
    mv "${out}.tmp" "${out}"
    rm "${raw}.txt"
}

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
        cat "${trial_files[@]}" | sort -k3 -gr | head -n ${K_TOP} > "${agg}.tmp"
        mv "${agg}.tmp" "${agg}"
    else
        echo "Skipping agg ${ds} min${t} trials${r} -- already done" >&2
    fi

    ## Evaluate against EXACT at each h.
    for h in "${H_VALUES[@]}"
    do
        local tsv="${OUT_DIR}/${ds}/eval/${ds}.slider++.agg.min${t}.trials${r}.k${K_TOP}.h${h}.tsv"
        if [ ! -f "${tsv}" ]
        then
            echo "Evaluating ${ds} min${t} trials${r} h${h}" >&2
            ${CMM_EVAL} "${fa}" "${int}" \
                "${agg}" "${exact}" \
                ${h} "${tsv}.tmp" /dev/null
            mv "${tsv}.tmp" "${tsv}"
        else
            echo "Skipping eval ${ds} min${t} trials${r} h${h} -- already done" >&2
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
            run_slider ${ds} ${trial} ${t}
        done
    done
    for t in "${TIMES[@]}"
    do
        agg_and_eval ${ds} ${t} ${TRIALS}
    done
done

## RQ2: vary number of trials (max time fixed at 5)
echo "Running experiment for RQ2..." >&2
TIME=5
TRIAL_COUNTS=(5 10 15 20 25 30)
MAX_TRIALS=30
for ds in "${DATASETS[@]}"
do
    for trial in $(seq 1 ${MAX_TRIALS})
    do
        run_slider ${ds} ${trial} ${TIME}
    done
    for r in "${TRIAL_COUNTS[@]}"
    do
        agg_and_eval ${ds} ${TIME} ${r}
    done
done

echo "All experiments done." >&2
