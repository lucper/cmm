#!/usr/bin/bash

set -euo pipefail

usage() {
    cat >&2 <<EOF
Usage: $0 -D <dataset> -f <fasta> -i <ints> -r <R> -t <T> -l <L> -d <D> -k <K> -K <K_TOP> -H <h> -c <C>

Options:
  -D <name>  dataset name
  -f <path>  FASTA file
  -i <path>  interactions file
  -r <int>   number of trials/runs
  -t <int>   max convergence time (min)
  -l <int>   motif length
  -d <int>   number of wildcards
  -k <int>   SLIDER output size K
  -K <int>   kept top pairs after aggregation
  -H <int>   similarity proximity h
  -c <int>   recall cutoff C in [0,100]
  -h         show this help and exit
EOF
}

DATASET="" FASTA="" INTS="" R="" T="" L="" D_WC="" K="" K_TOP="" H="" C=""

while getopts "D:f:i:r:t:l:d:k:K:H:c:h" opt
do
    case ${opt} in
        D) DATASET=${OPTARG} ;;
        f) FASTA=${OPTARG} ;;
        i) INTS=${OPTARG} ;;
        r) R=${OPTARG} ;;
        t) T=${OPTARG} ;;
        l) L=${OPTARG} ;;
        d) D_WC=${OPTARG} ;;
        k) K=${OPTARG} ;;
        K) K_TOP=${OPTARG} ;;
        H) H=${OPTARG} ;;
        c) C=${OPTARG} ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

declare -A FLAG=(
    [DATASET]=-D [FASTA]=-f [INTS]=-i [R]=-r [T]=-t [L]=-l
    [D_WC]=-d [K]=-k [K_TOP]=-K [H]=-H [C]=-c
)
missing=""
for var in DATASET FASTA INTS R T L D_WC K K_TOP H C
do
    [ -z "${!var}" ] && missing="${missing} ${FLAG[$var]}"
done
if [ -n "${missing}" ]
then
    echo "Error: missing required option(s):${missing}" >&2
    usage
    exit 1
fi
for var in R T L D_WC K K_TOP H C
do
    if ! [[ "${!var}" =~ ^[0-9]+$ ]]
    then
        echo "Error: ${FLAG[$var]} must be a non-negative integer, got '${!var}'." >&2
        exit 1
    fi
done

FASTA=$(realpath "${FASTA}")
INTS=$(realpath "${INTS}")

echo "=== Running experiment with parameters: DATASET=${DATASET} R=${R} T=${T} L=${L} D=${D_WC} K=${K} K_TOP=${K_TOP} H=${H} C=${C} ===" >&2

echo "=== Stage 010-soln: SLIDER runs ===" >&2
for i in $(seq 1 ${R}); do
    make -C 010-soln all \
        DATASET=${DATASET} FASTA=${FASTA} INTS=${INTS} I=${i} T=${T} L=${L} D=${D_WC} K=${K}
done

echo "=== Stage 020-agg: Aggregation ===" >&2
make -C 020-agg all \
    DATASET=${DATASET} R=${R} T=${T} L=${L} D=${D_WC} K=${K} K_TOP=${K_TOP}

echo "=== Stage 025-dedup: Deduplication ===" >&2
make -C 025-dedup slider \
    DATASET=${DATASET} FASTA=${FASTA} INTS=${INTS} R=${R} T=${T} L=${L} D=${D_WC} K_TOP=${K_TOP}
make -C 025-dedup cmm \
    DATASET=${DATASET} FASTA=${FASTA} INTS=${INTS} L=${L} D=${D_WC} K_TOP=${K_TOP}

echo "=== Stage 030-eval: Similarity evaluation ===" >&2
make -C 030-eval all \
    DATASET=${DATASET} FASTA=${FASTA} INTS=${INTS} R=${R} T=${T} L=${L} D=${D_WC} K_TOP=${K_TOP} H=${H}

echo "=== Stage 035-eval-recall: Recall against exact ===" >&2
make -C 035-eval-recall all \
    DATASET=${DATASET} R=${R} T=${T} L=${L} D=${D_WC} K_TOP=${K_TOP} H=${H} C=${C}
