#!/usr/bin/bash

set -euo pipefail

usage() {
    cat >&2 <<EOF
Usage: $0 -D <dataset> -f <fasta> -i <ints> -r <R> -t <T> -l <L> -d <Ds> -k <K> -K <K_TOP> -H <h> -c <C> -p <P>

Options:
  -D <name>  dataset name
  -f <path>  FASTA file
  -i <path>  interactions file
  -r <int>   number of trials/runs
  -t <int>   max convergence time (min)
  -l <int>   motif length
  -d <list>  numbers of wildcards (e.g. "3 5"); the stages run once per value
  -k <int>   SLIDER output size K
  -K <int>   kept top pairs after aggregation
  -H <int>   similarity proximity h
  -c <int>   recall cutoff C in [0,100]
  -p <int>   threads for the exact miner (affects speed only)
  -h         show this help and exit
EOF
}

DATASET="" FASTA="" INTS="" R="" T="" L="" DS="" K="" K_TOP="" H="" C="" CMM_THREADS=""

while getopts "D:f:i:r:t:l:d:k:K:H:c:p:h" opt
do
    case ${opt} in
        D) DATASET=${OPTARG} ;;
        f) FASTA=${OPTARG} ;;
        i) INTS=${OPTARG} ;;
        r) R=${OPTARG} ;;
        t) T=${OPTARG} ;;
        l) L=${OPTARG} ;;
        d) DS=${OPTARG} ;;
        k) K=${OPTARG} ;;
        K) K_TOP=${OPTARG} ;;
        H) H=${OPTARG} ;;
        c) C=${OPTARG} ;;
        p) CMM_THREADS=${OPTARG} ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

declare -A FLAG=(
    [DATASET]=-D [FASTA]=-f [INTS]=-i [R]=-r [T]=-t [L]=-l
    [DS]=-d [K]=-k [K_TOP]=-K [H]=-H [C]=-c [CMM_THREADS]=-p
)
missing=""
for var in DATASET FASTA INTS R T L DS K K_TOP H C CMM_THREADS
do
    [ -z "${!var}" ] && missing="${missing} ${FLAG[$var]}"
done
if [ -n "${missing}" ]
then
    echo "Error: missing required option(s):${missing}" >&2
    usage
    exit 1
fi
for var in R T L K K_TOP H C CMM_THREADS
do
    if ! [[ "${!var}" =~ ^[0-9]+$ ]]
    then
        echo "Error: ${FLAG[$var]} must be a non-negative integer, got '${!var}'." >&2
        exit 1
    fi
done
for d in ${DS}
do
    if ! [[ "${d}" =~ ^[0-9]+$ ]] || [ "${d}" -ge "${L}" ]
    then
        echo "Error: -d must be a list of integers in [0,${L}), got '${DS}'." >&2
        exit 1
    fi
done

FASTA=$(realpath "${FASTA}")
INTS=$(realpath "${INTS}")

for D_WC in ${DS}
do
    echo "=== Running experiment with parameters: DATASET=${DATASET} R=${R} T=${T} L=${L} D=${D_WC} K=${K} K_TOP=${K_TOP} H=${H} C=${C} ===" >&2

    echo "=== Stage 010-soln-slider: SLIDER runs ===" >&2
    for i in $(seq 1 ${R}); do
        make -C 010-soln-slider all \
            DATASET=${DATASET} FASTA=${FASTA} INTS=${INTS} I=${i} T=${T} L=${L} D=${D_WC} K=${K}
    done

    echo "=== Stage 010-soln-cmm: exact solution ===" >&2
    make -C 010-soln-cmm all \
        DATASET=${DATASET} FASTA=${FASTA} INTS=${INTS} L=${L} D=${D_WC} K_TOP=${K_TOP} THREADS=${CMM_THREADS}

    echo "=== Stage 020-agg: Aggregation ===" >&2
    make -C 020-agg all \
        DATASET=${DATASET} R=${R} T=${T} L=${L} D=${D_WC} K=${K} K_TOP=${K_TOP}

    echo "=== Stage 025-dedup: Deduplication ===" >&2
    make -C 025-dedup slider \
        DATASET=${DATASET} FASTA=${FASTA} INTS=${INTS} R=${R} T=${T} L=${L} D=${D_WC} K_TOP=${K_TOP}
    make -C 025-dedup cmm \
        DATASET=${DATASET} FASTA=${FASTA} INTS=${INTS} L=${L} D=${D_WC} K_TOP=${K_TOP}

    echo "=== Stage 030-score-curve: Score curve ===" >&2
    make -C 030-score-curve all \
        DATASET=${DATASET} R=${R} T=${T} L=${L} D=${D_WC} K_TOP=${K_TOP}

    echo "=== Stage 030-eval: Similarity evaluation ===" >&2
    make -C 030-eval all \
        DATASET=${DATASET} FASTA=${FASTA} INTS=${INTS} R=${R} T=${T} L=${L} D=${D_WC} K_TOP=${K_TOP} H=${H}

    echo "=== Stage 035-eval-density: Similarity density plot ===" >&2
    make -C 035-eval-density all \
        DATASET=${DATASET} R=${R} T=${T} L=${L} D=${D_WC} K_TOP=${K_TOP} H=${H}

    echo "=== Stage 035-eval-coverage: Coverage of exact by heuristics ===" >&2
    make -C 035-eval-coverage all \
        DATASET=${DATASET} R=${R} T=${T} L=${L} D=${D_WC} K_TOP=${K_TOP} H=${H} C=${C}
done
