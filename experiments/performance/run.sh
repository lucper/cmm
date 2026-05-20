#!/usr/bin/bash

set -euo pipefail

DATA="../data/aj-human/BioPlex_293T.seq"
GENERATE="../generate_instance.py"
CMM_PERF="./cmm_perf"
VENV="venv"
SEED=42

OUT_RQ1="rq1.tsv"
OUT_RQ2="rq2.tsv"
OUT_RQ3="rq3.tsv"

MAX_THREADS=64
TRIALS=5

## Setup virtual environment
if [ ! -d "${VENV}" ]; then
    echo "Creating virtual environment..." >&2
    python3 -m venv ${VENV}
    ${VENV}/bin/pip install -q -r requirements.txt
fi
source ${VENV}/bin/activate

## Start plot watcher
#echo "Starting plot watcher..." >&2
#nohup watch -n 300 make -C plots > plots/watch.log 2>&1 &
#WATCH_PID=$!
#trap "kill ${WATCH_PID}" EXIT

## Generate artificial instances from real data
echo "Generating instances..." >&2
for ed in 05 10 15 20 25 30
do
    mkdir -p instances/ed${ed}
    for v in 100 200 400 800 1600
    do
        fa="instances/ed${ed}/sampled_human_V${v}_ed${ed}.fa"
        int="instances/ed${ed}/sampled_human_V${v}_ed${ed}.int"
        if [ ! -f "${fa}" ] || [ ! -f "${int}" ]; then
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
    local ed=$7
    local threads=$8

    local key=$(basename ${fa} .fa)

    if [ -f "${tsv}" ] && grep -q "^${trial}	${key}" ${tsv}; then
        echo "Skipping trial=${trial} ${key} ell=${ell} d=${d} threads=${threads} -- already done" >&2
        return
    fi

    echo "Running trial=${trial} ${key} ell=${ell} d=${d} threads=${threads}" >&2
    ${CMM_PERF} ${fa} ${int} ${ell} ${d} 0.${ed} ${threads} \
        | awk -v trial=${trial} '{print trial "\t" $0}' >> ${tsv}
}

## RQ1
echo "Running experiment for RQ1..." >&2
ED=05
for trial in $(seq 1 ${TRIALS})
do
    for v in 100 200 400 800 1600
    do
        fa="instances/ed${ED}/sampled_human_V${v}_ed${ED}.fa"
        int="instances/ed${ED}/sampled_human_V${v}_ed${ED}.int"

        ell=5
        for d in 0 1 2
        do
            run_cmm_perf ${OUT_RQ1} ${trial} ${fa} ${int} ${ell} ${d} ${ED} ${MAX_THREADS}
        done

        ell=8
        for d in 0 1 2 3 4
        do
            run_cmm_perf ${OUT_RQ1} ${trial} ${fa} ${int} ${ell} ${d} ${ED} ${MAX_THREADS}
        done
    done
done

## RQ2
echo "Running experiment for RQ2..." >&2
V=100
for trial in $(seq 1 ${TRIALS})
do
    for ed in 05 10 15 20 25 30
    do
        fa="instances/ed${ed}/sampled_human_V${V}_ed${ed}.fa"
        int="instances/ed${ed}/sampled_human_V${V}_ed${ed}.int"

        ell=5
        for d in 0 1 2
        do
            run_cmm_perf ${OUT_RQ2} ${trial} ${fa} ${int} ${ell} ${d} ${ed} ${MAX_THREADS}
        done

        ell=8
        for d in 0 1 2 3 4
        do
            run_cmm_perf ${OUT_RQ2} ${trial} ${fa} ${int} ${ell} ${d} ${ed} ${MAX_THREADS}
        done
    done
done

## RQ3
echo "Running experiment for RQ3..." >&2
ED=05
V=100
ELL=8
D=4
for trial in $(seq 1 ${TRIALS})
do
    for t in 1 2 4 8 16 32 64
    do
        fa="instances/ed${ED}/sampled_human_V${V}_ed${ED}.fa"
        int="instances/ed${ED}/sampled_human_V${V}_ed${ED}.int"
        run_cmm_perf ${OUT_RQ3} ${trial} ${fa} ${int} ${ELL} ${D} ${ED} ${t}
    done
done

echo "All experiments done." >&2
