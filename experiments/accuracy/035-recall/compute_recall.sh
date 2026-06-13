#!/usr/bin/bash

set -euo pipefail

if [ "$#" -ne 3 ]
then
    echo "Usage: $0 <cmm_eval_h0_dedup_tsv> <cmm_dedup_file> <out>" >&2
    exit 1
fi

eval_tsv=$1
cmm_dedup=$2
out=$3

baseline_set=$(wc -l < "${cmm_dedup}")
competitor_set=$(wc -l < "${eval_tsv}")

if [ "${baseline_set}" -eq 0 ]
then
    echo "Error: baseline set is empty (cmm_dedup is empty: ${cmm_dedup})" >&2
    exit 1
fi

matches=$(awk -v k="${baseline_set}" '$1 <= k && $8 + 0 == 1' "${eval_tsv}" | wc -l)
competitor_set_trunc=$(( competitor_set < baseline_set ? competitor_set : baseline_set ))
recall=$(awk -v m="${matches}" -v k="${baseline_set}" 'BEGIN{printf "%.4f", m/k}')

printf "%d\t%d\t%d\t%s\n" "${baseline_set}" "${competitor_set_trunc}" "${matches}" "${recall}" > "${out}"
