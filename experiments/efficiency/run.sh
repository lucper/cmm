#!/usr/bin/bash

set -euo pipefail

DATA="../string_data_cleaned/human.s700.cleaned.fa"
GENERATE="../generate_instance.py"
CMM_PERF="./cmm_perf"
VENV="venv"
SEED=42

OUT_RQ1="rq1_out.tsv"
OUT_RQ2="rq2_out.tsv"
OUT_RQ3="rq3_out.tsv"

MAX_THREADS=64
TRIALS=5
K=1000

## Setup virtual environment
if [ ! -d "${VENV}" ]
then
    echo "Creating virtual environment..." >&2
    python3 -m venv ${VENV}
    ${VENV}/bin/pip install -q -r requirements.txt
fi
source ${VENV}/bin/activate

## Generate artificial instances from real data
echo "Generating instances..." >&2
for ed in 05 10 15 20 25 30
do
    mkdir -p instances/ed${ed}
    for v in 50 100 200 400 800 1600 3200
    do
        fa="instances/ed${ed}/sampled_human_V${v}_ed${ed}.fa"
        int="instances/ed${ed}/sampled_human_V${v}_ed${ed}.int"
        if [ ! -f "${fa}" ] || [ ! -f "${int}" ]
        then
            ${GENERATE} \
                --n ${v} --density 0.${ed} --seed ${SEED} \
                --out-nodes ${fa} \
                --out-edges ${int} \
                --fasta ${DATA}
        fi
    done
done

run_cmm_perf() {
    local tsv=$1
    local trial=$2
    local fa=$3
    local int=$4
    local ell=$5
    local d=$6
    local threads=$7

    local key=$(basename ${fa} .fa)

    # TSV columns after awk-prepended trial:
    # 1:trial 2:instance 3:edge_density 4:V 5:N 6:ell 7:d
    # 8:max_rank 9:min_rank 10:avg_rank 11:pruning_cnt
    # 12:num_threads_requested 13:num_threads_spawned 14:time_ms 15:peak_ram_kb
    local uniq_pattern="${trial}|${key}|${ell}|${d}|${threads}"
    if [ -f "${tsv}" ] && (cut -f 1,2,6,7,12 ${tsv} | tr '\t' '|' | grep -qF "${uniq_pattern}")
    then
        echo "Skipping trial=${trial} ${key} ell=${ell} d=${d} threads=${threads} -- already done" >&2
        return
    fi

    echo "Running trial=${trial} ${key} ell=${ell} d=${d} threads=${threads}" >&2
    ${CMM_PERF} ${fa} ${int} ${ell} ${d} ${K} ${threads} /dev/null \
        | awk -v trial=${trial} '{print trial "\t" $0}' >> ${tsv}
}

## RQ1
echo "Running experiment for RQ1..." >&2
ED=05
for trial in $(seq 1 ${TRIALS})
do
    for v in 100 200 400 800 1600 3200
    do
        fa="instances/ed${ED}/sampled_human_V${v}_ed${ED}.fa"
        int="instances/ed${ED}/sampled_human_V${v}_ed${ED}.int"

        ell=5
        for d in 0 1 2 3 4
        do
            run_cmm_perf ${OUT_RQ1} ${trial} ${fa} ${int} ${ell} ${d} ${MAX_THREADS}
        done

        ell=8
        for d in 0 1 2 3 4 5 6 7
        do
            run_cmm_perf ${OUT_RQ1} ${trial} ${fa} ${int} ${ell} ${d} ${MAX_THREADS}
        done
    done
done

## RQ2
echo "Running experiment for RQ2..." >&2
V=3200
for trial in $(seq 1 ${TRIALS})
do
    for ed in 05 10 15 20 25 30
    do
        fa="instances/ed${ed}/sampled_human_V${V}_ed${ed}.fa"
        int="instances/ed${ed}/sampled_human_V${V}_ed${ed}.int"

        ell=5
        for d in 0 1 2 3 4
        do
            run_cmm_perf ${OUT_RQ2} ${trial} ${fa} ${int} ${ell} ${d} ${MAX_THREADS}
        done

        ell=8
        for d in 0 1 2 3 4 5 6 7
        do
            run_cmm_perf ${OUT_RQ2} ${trial} ${fa} ${int} ${ell} ${d} ${MAX_THREADS}
        done
    done
done

## RQ3
echo "Running experiment for RQ3..." >&2
ED=05
V=3200
ELL=5
D=3
for trial in $(seq 1 ${TRIALS})
do
    for t in 1 2 4 8 16 32 64
    do
        fa="instances/ed${ED}/sampled_human_V${V}_ed${ED}.fa"
        int="instances/ed${ED}/sampled_human_V${V}_ed${ED}.int"
        run_cmm_perf ${OUT_RQ3} ${trial} ${fa} ${int} ${ELL} ${D} ${t}
    done
done

echo "All experiments done." >&2
