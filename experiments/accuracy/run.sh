#!/usr/bin/bash

set -euo pipefail

DATASETS=(ecoli yeast fly plant) # add human
TRIALS=5
MAXTIMES=(5 10 15 20 25 30)
ELL=8
D=3
K_SLIDER=10000
K_TOP=1000
H_VALUES=(0 ${ELL})

SEQ_DIR="../string_data_cleaned" # change this
CMM_DIR="cmm_out"                # change this
OUT_DIR="slider_out"             # change this

CMM_EVAL=./cmm_eval
CMM_DEDUP=./cmm_dedup
SLIDER_JAR="../competitors/SliderLight/dist/SliderLight.jar"

## sanity checks: fail if anything's missing
for tool in "${CMM_EVAL}" "${CMM_DEDUP}"
do
    [ -x "${tool}" ] || { echo "Error: ${tool} not found or not executable" >&2; exit 1; }
done
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

## slider output
for ds in "${DATASETS[@]}"
do
    mkdir -p "${OUT_DIR}/${ds}" "${OUT_DIR}/${ds}/agg" "${OUT_DIR}/${ds}/eval"
done

echo "=== Stage 1: SliderLight runs ===" >&2
for ds in "${DATASETS[@]}"
do
    fa="${SEQ_DIR}/${ds}.s700.cleaned.fa"
    int="${SEQ_DIR}/${ds}.s700.cleaned.int"
    for i in $(seq 1 ${TRIALS})
    do
        for t in "${MAXTIMES[@]}"
        do
            out="${OUT_DIR}/${ds}/${ds}.slider++.trial${i}.min${t}.k${K_SLIDER}.out"
            if [ -f "${out}" ]
            then
                echo "Skip:    SLIDER ${ds} trial${i} min${t}" >&2
                continue
            fi

            echo "Run:     SLIDER ${ds} trial${i} min${t} (up to ${t} min)" >&2
            raw="${out}.raw"
            java -cp "${SLIDER_JAR}" Framework.Framework \
                -l ${ELL} -d ${D} \
                -seq "${fa}" \
                -int "${int}" \
                -o "${raw}" \
                -m slider++ -a ${K_SLIDER} -st x2 -min ${t}
            tr '-' ' ' < "${raw}.txt" > "${out}.tmp"
            mv "${out}.tmp" "${out}"
            rm "${raw}.txt"
        done
    done
done

## concat and sort outputs from multiple trials
echo "=== Stage 2: Aggregations ===" >&2
for ds in "${DATASETS[@]}"
do
    for t in "${MAXTIMES[@]}"
    do
        trial_files=()
        for i in $(seq 1 ${TRIALS})
        do
            trial_files+=("${OUT_DIR}/${ds}/${ds}.slider++.trial${i}.min${t}.k${K_SLIDER}.out")
        done

        ## No-dedup
        out_all="${OUT_DIR}/${ds}/agg/${ds}.slider++.agg.min${t}.all.out"
        if [ ! -f "${out_all}" ]
        then
            echo "Agg:     ${ds} min${t} sorted aggregate" >&2
            cat "${trial_files[@]}" | sort -k3 -gr > "${out_all}.tmp"
            mv "${out_all}.tmp" "${out_all}"
        else
            echo "Skip:    Agg ${ds} min${t} full" >&2
        fi
    done
done

echo "=== Stage 3: Dedup + truncate ===" >&2
for ds in "${DATASETS[@]}"
do
    fa="${SEQ_DIR}/${ds}.s700.cleaned.fa"
    int="${SEQ_DIR}/${ds}.s700.cleaned.int"
    for t in "${MAXTIMES[@]}"
    do
        all_out="${OUT_DIR}/${ds}/agg/${ds}.slider++.agg.min${t}.all.out"
        dedup_full="${OUT_DIR}/${ds}/agg/${ds}.slider++.agg.min${t}.all.dedup.out"

        if [ ! -f "${dedup_full}" ]; then
            echo "Dedup:   ${ds} min${t} (full aggregate)" >&2
            ${CMM_DEDUP} "${fa}" "${int}" "${all_out}" 0 "${dedup_full}.tmp"
            mv "${dedup_full}.tmp" "${dedup_full}"
        else
            echo "Skip:    Dedup ${ds} min${t}" >&2
        fi
    done
done

echo "=== Stage 4: cmm_eval ===" >&2
for ds in "${DATASETS[@]}"
do
    fa="${SEQ_DIR}/${ds}.s700.cleaned.fa"
    int="${SEQ_DIR}/${ds}.s700.cleaned.int"
    exact="${CMM_DIR}/${ds}.cmm.out"

    for t in "${MAXTIMES[@]}"
    do
        in_nd="${OUT_DIR}/${ds}/agg/${ds}.slider++.agg.min${t}.k${K_TOP}.out"
        if [ ! -f "${in_nd}" ]
        then
            head -n${K_TOP} "${OUT_DIR}/${ds}/agg/${ds}.slider++.agg.min${t}.all.out"  > "${in_nd}"
        fi

        in_d="${OUT_DIR}/${ds}/agg/${ds}.slider++.agg.min${t}.k${K_TOP}.dedup.out"
        if [ ! -f "${in_d}" ]
        then
            ## note that cmm_dedup output has a header, so we take head -n K_TOP+1
            head -n$((K_TOP + 1)) "${OUT_DIR}/${ds}/agg/${ds}.slider++.agg.min${t}.all.dedup.out" > "${in_d}"
        fi

        for h in "${H_VALUES[@]}"
        do
            # no-dedup evaluation
            tsv_nd="${OUT_DIR}/${ds}/eval/${ds}.slider++.agg.min${t}.k${K_TOP}.h${h}.tsv"
            if [ ! -f "${tsv_nd}" ]
            then
                echo "Eval:    ${ds} min${t} h${h} no-dedup" >&2
                ${CMM_EVAL} "${fa}" "${int}" \
                    "${in_nd}" "${exact}" \
                    ${h} "${tsv_nd}.tmp" /dev/null
                mv "${tsv_nd}.tmp" "${tsv_nd}"
            else
                echo "Skip:    Eval ${ds} min${t} h${h} no-dedup" >&2
            fi

            # dedup evaluation
            tsv_d="${OUT_DIR}/${ds}/eval/${ds}.slider++.agg.min${t}.k${K_TOP}.h${h}.dedup.tsv"
            if [ ! -f "${tsv_d}" ]
            then
                echo "Eval:    ${ds} min${t} h${h} dedup" >&2
                ${CMM_EVAL} "${fa}" "${int}" \
                    "${in_d}" "${exact}" \
                    ${h} "${tsv_d}.tmp" /dev/null
                mv "${tsv_d.tmp}" "${tsv_d}"
            else
                echo "Skip:    Eval ${ds} min${t} h${h} dedup" >&2
            fi
        done
    done
done

echo "=== All done ===" >&2
