#!/usr/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <main_algo_output.tsv>"
    exit 1
fi

cat $1 | awk 'NR > 1 && $1 < $2 { tmp=$1; $1=$2; $2=tmp } NR > 1 { print $1, $2, $3 }' OFS='\t' | sort
