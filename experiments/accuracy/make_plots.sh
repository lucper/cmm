#!/usr/bin/bash

set -euo pipefail

declare -A DATASETS=( ["ecoli"]="E. coli" ["yeast"]="S. cerevisiae" ["fly"]="D. melanogaster" ["plant"]="A. thaliana" )

ELL=8
H_VALUES=(0 ${ELL})
K_TOP=1000
THRESHOLD=0.5

OUT_DIR="slider_out"
PLOT_OUT_DIR="plots"
PLOT_SCRIPT="./plot_top_n_recall.py"


mkdir -p "${PLOT_OUT_DIR}"

[ -f "${PLOT_SCRIPT}" ] || { echo "Error: ${PLOT_SCRIPT} not found" >&2; exit 1; }

## Build INPUTS array for a given (rq, h, suffix) cell.
## Sets global INPUTS; sets COUNT_FOUND and COUNT_MISSING for diagnostics.
build_inputs() {
    local rq=$1
    local h=$2
    local suffix=$3  # "" / ".dedup" / ".dedupall"

    INPUTS=()
    COUNT_FOUND=0
    COUNT_MISSING=0

    local ds t r path
    for ds in "${!DATASETS[@]}"
    do
        if [ "${rq}" = "rq1" ]
        then
            for t in 5 10 15 20 25 30
            do
                path="${OUT_DIR}/${ds}/eval/${ds}.slider++.agg.min${t}.trials5.k${K_TOP}.h${h}${suffix}.tsv"
                if [ -f "${path}" ]
                then
                    INPUTS+=("${DATASETS[${ds}]}:${t}:${path}")
                    COUNT_FOUND=$((COUNT_FOUND + 1))
                else
                    COUNT_MISSING=$((COUNT_MISSING + 1))
                fi
            done
        else  # rq2
            for r in 5 10 15 20 25 30
            do
                path="${OUT_DIR}/${ds}/eval/${ds}.slider++.agg.min5.trials${r}.k${K_TOP}.h${h}${suffix}.tsv"
                if [ -f "${path}" ]
                then
                    INPUTS+=("${DATASETS[${ds}]}:${r}:${path}")
                    COUNT_FOUND=$((COUNT_FOUND + 1))
                else
                    COUNT_MISSING=$((COUNT_MISSING + 1))
                fi
            done
        fi
    done
}

## Loop over (RQ, h, dedup-mode) and produce each plot if input exists.
for rq in rq1 rq2
do
    if [ "${rq}" = "rq1" ]
    then
        xlabel="Max convergence time (min)"
    else
        xlabel="Number of trials"
    fi

    for h in "${H_VALUES[@]}"
    do
        ## Three dedup variants encoded as "tag:filename_suffix" pairs.
        for variant in "nodedup:" "headdedup:.dedup" "dedupall:.dedupall"
        do
            tag="${variant%:*}"
            suffix="${variant#*:}"

            build_inputs "${rq}" "${h}" "${suffix}"

            if [ ${COUNT_FOUND} -eq 0 ]
            then
                echo "Skip:    ${rq} h${h} ${tag} -- 0/${COUNT_MISSING} input TSVs present" >&2
                continue
            fi

            out="${PLOT_OUT_DIR}/acc_${rq}_h${h}_${tag}.pdf"
            echo "Plot:    ${out} (${COUNT_FOUND} input TSVs, ${COUNT_MISSING} missing)" >&2

            ./venv/bin/python3 "${PLOT_SCRIPT}" \
                --input "${INPUTS[@]}" \
                --output "${out}" \
                --x-label "${xlabel}" \
                --top-n ${K_TOP} \
                --threshold ${THRESHOLD}
        done
    done
done

echo "=== Plots in ${PLOT_OUT_DIR}/ ===" >&2
ls -1 "${PLOT_OUT_DIR}"
