#!/usr/bin/bash

cat $1 | awk 'NR > 1 && $1 < $2 { tmp=$1; $1=$2; $2=tmp } NR > 1 { print $1, $2, $3 }' OFS='\t' | sort
