#!/usr/bin/bash
#
# 1. fetch_string.sh   downloads STRING sequences + links at a score cutoff
# 2. clean_dataset.py  drops non-amino-acid records, dangling edges, isolates
# 3. generate_instance.py  samples artificial instances at target densities
#
# Usage:
#   ./run.sh

set -euo pipefail

cd "$(dirname "$0")"

## Parameters
THRESHOLD=700
TAXID=9606
ORGANISM=Homo_sapiens
SEED=42
DENSITIES="05 10 15 20 25 30"
NODE_COUNTS="50 100 200 400 800 1600 3200"

FETCH="./scripts/fetch_string.sh"
CLEAN="./scripts/clean_dataset.py"
GENERATE="./scripts/generate_instance.py"

STRING_DIR="string_data/${TAXID}_${ORGANISM}"
RAW_FA="${STRING_DIR}/${TAXID}.sequences.s${THRESHOLD}.fa"
RAW_INT="${STRING_DIR}/${TAXID}.links.s${THRESHOLD}.int"

CLEAN_DIR="real_instances"
CLEAN_FA="${CLEAN_DIR}/${TAXID}.sequences.s${THRESHOLD}.clean.fa"
CLEAN_INT="${CLEAN_DIR}/${TAXID}.links.s${THRESHOLD}.clean.int"

INSTANCE_DIR="artificial_instances"

## 1. Fetch raw data from STRING
echo "Fetching STRING data (cutoff ${THRESHOLD})..." >&2
if [ ! -s "${RAW_FA}" ] || [ ! -s "${RAW_INT}" ]
then
    ${FETCH} ${THRESHOLD}
else
    echo "  ${RAW_FA} already present, skipping" >&2
fi

## 2. Clean: drop dubious sequences, dangling edges, isolated nodes
echo "Cleaning dataset..." >&2
mkdir -p ${CLEAN_DIR}
if [ ! -s "${CLEAN_FA}" ] || [ ! -s "${CLEAN_INT}" ]
then
    ${CLEAN} \
        --fasta ${RAW_FA} \
        --edges ${RAW_INT} \
        --out-nodes ${CLEAN_FA} \
        --out-edges ${CLEAN_INT}
else
    echo "  ${CLEAN_FA} already present, skipping" >&2
fi

## 3. Generate artificial instances from cleaned data
echo "Generating instances..." >&2
for ed in ${DENSITIES}
do
    mkdir -p ${INSTANCE_DIR}/ed${ed}
    for v in ${NODE_COUNTS}
    do
        fa="${INSTANCE_DIR}/ed${ed}/sampled_human_V${v}_ed${ed}.fa"
        int="${INSTANCE_DIR}/ed${ed}/sampled_human_V${v}_ed${ed}.int"
        if [ ! -f "${fa}" ] || [ ! -f "${int}" ]
        then
            echo "  V=${v} ed=0.${ed}" >&2
            ${GENERATE} \
                --fasta ${CLEAN_FA} \
                --n ${v} --density 0.${ed} --seed ${SEED} \
                --out-nodes ${fa} \
                --out-edges ${int}
        fi
    done
done

echo "Done. Instances written to ${INSTANCE_DIR}/ed<NN>/" >&2
