#!/usr/bin/bash
#
# Layout produced under this directory:
#   string_data/<taxid>_<name>/   raw + per-cutoff .fa/.int (from fetch_string.sh)
#   real_instances/               cleaned .fa/.int, one pair per organism
#   artificial_instances/ed<NN>/  sampled instances
#
# Usage:
#   ./run.sh

set -euo pipefail

cd "$(dirname "$0")"

## Parameters
THRESHOLD=700
SEED=42

# Organisms to fetch and clean, as "taxid:name" (passed through to fetch_string.sh).
# Instances are generated from INSTANCE_TAXID only.
ORGANISMS=(
    "9606:Homo_sapiens"
    "3702:Arabidopsis_thaliana"
    "4932:Saccharomyces_cerevisiae"
    "511145:Escherichia_coli_K12_MG1655"
    "7227:Drosophila_melanogaster"
)

source ../params.conf

SCRIPT_DIR="scripts"
FETCH="${SCRIPT_DIR}/fetch_string.sh"
CLEAN="${SCRIPT_DIR}/clean_dataset.py"
GENERATE="${SCRIPT_DIR}/generate_instance.py"

STRING_DIR="string_data"
CLEAN_DIR="real_instances"
INSTANCE_DIR="artificial_instances"

# Path helpers, so the layout is defined in exactly one place.
raw_fa()    { echo "${STRING_DIR}/$1_$2/$1.sequences.s${THRESHOLD}.fa"; }
raw_int()   { echo "${STRING_DIR}/$1_$2/$1.links.s${THRESHOLD}.int"; }
clean_fa()  { echo "${CLEAN_DIR}/$1.sequences.s${THRESHOLD}.clean.fa"; }
clean_int() { echo "${CLEAN_DIR}/$1.links.s${THRESHOLD}.clean.int"; }

## 1. Fetch raw data from STRING (all organisms)
echo "Fetching STRING data (cutoff ${THRESHOLD})..." >&2
need_fetch=0
for entry in "${ORGANISMS[@]}"
do
    taxid="${entry%%:*}"
    name="${entry#*:}"
    if [ ! -s "$(raw_fa ${taxid} ${name})" ] || [ ! -s "$(raw_int ${taxid} ${name})" ]
    then
        need_fetch=1
        break
    fi
done

if [ ${need_fetch} -eq 1 ]
then
    OUTDIR="${STRING_DIR}" ${FETCH} ${THRESHOLD} "${ORGANISMS[@]}"
else
    echo "  all organisms already fetched, skipping" >&2
fi

## 2. Clean each organism: drop dubious sequences, dangling edges, isolated nodes
echo "Cleaning datasets..." >&2
mkdir -p ${CLEAN_DIR}
for entry in "${ORGANISMS[@]}"
do
    taxid="${entry%%:*}"
    name="${entry#*:}"

    fa_in=$(raw_fa ${taxid} ${name})
    int_in=$(raw_int ${taxid} ${name})
    fa_out=$(clean_fa ${taxid})
    int_out=$(clean_int ${taxid})

    if [ -s "${fa_out}" ] && [ -s "${int_out}" ]
    then
        echo "  ${name}: already cleaned, skipping" >&2
        continue
    fi

    echo "  ${name} (taxid ${taxid}):" >&2
    ${CLEAN} \
        --fasta ${fa_in} \
        --edges ${int_in} \
        --out-nodes ${fa_out} \
        --out-edges ${int_out}
done

## 3. Generate artificial instances (Homo sapiens only)
echo "Generating instances from taxid ${INSTANCE_TAXID}..." >&2
INSTANCE_SOURCE=$(clean_fa ${INSTANCE_TAXID})

num_seqs=$(grep -c '^>' "${INSTANCE_SOURCE}")
max_v=$(echo ${NODE_COUNTS} | tr ' ' '\n' | sort -n | tail -1)
if [ "${num_seqs}" -lt "${max_v}" ]
then
    echo "Error: ${INSTANCE_SOURCE} has ${num_seqs} sequences, but NODE_COUNTS needs ${max_v}." >&2
    exit 1
fi

for ed in ${DENSITIES}
do
    mkdir -p ${INSTANCE_DIR}/ed${ed}
    for v in ${NODE_COUNTS}
    do
        fa="${INSTANCE_DIR}/ed${ed}/sampled_${INSTANCE_TAXID}_V${v}_ed${ed}.fa"
        int="${INSTANCE_DIR}/ed${ed}/sampled_${INSTANCE_TAXID}_V${v}_ed${ed}.int"
        if [ ! -f "${fa}" ] || [ ! -f "${int}" ]
        then
            echo "  V=${v} ed=0.${ed}" >&2
            ${GENERATE} \
                --fasta ${INSTANCE_SOURCE} \
                --n ${v} --density 0.${ed} --seed ${SEED} \
                --out-nodes ${fa} \
                --out-edges ${int}
        fi
    done
done

echo "Done." >&2
echo "  cleaned data (all organisms): ${CLEAN_DIR}/" >&2
echo "  instances (${INSTANCE_TAXID}):  ${INSTANCE_DIR}/ed<NN>/" >&2
