#!/usr/bin/bash

DATA="../data/aj-human/BioPlex_293T.seq"

OUT_RQ1="rq1.tsv"
OUT_RQ2="rq2.tsv"
OUT_RQ3="rq3.tsv"

MAX_THREADS=64

## Generate artificial instances from real data
for ed in 05 10 15 20 25 30
do
    mkdir -p instances/ed${ed}
    for v in 100 200 400 800 1600
    do
        ../generate_instance.py \
            --n ${v} --density 0.${ed} --seed 42 \
            --out-nodes instances/ed${ed}/sampled_human_V${v}_ed${ed}.fa \
            --out-edges instances/ed${ed}/sampled_human_V${v}_ed${ed}.int \
            --fasta ${DATA}
    done
done

## RQ1: impact of V and N
#ED=05
#for V in 100 200 400 800 1600
#do
#    ell=5
#    for d in 0 1 2
#    do
#        ./cmm_perf \
#            instances/ed${ED}/sampled_human_V${V}_ed${ED}.fa \
#            instances/ed${ED}/sampled_human_V${V}_ed${ED}.int \
#            ${ell} ${d} 0.${ED} ${MAX_THREADS} >> ${OUT_RQ1}
#    done
#
#    ell=8
#    for d in 0 1 2 3 4
#    do
#        ./cmm_perf \
#            instances/ed{$ED}/sampled_human_V${V}_ed${ED}.fa \
#            instances/ed${ED}/sampled_human_V${V}_ed${ED}.int \
#            ${ell} ${d} 0.${ED} ${MAX_THREADS} >> ${OUT_RQ1}
#    done
#done

## RQ2: impact of edge density
V=100
for trial in 1 2 3 4 5
do
    for e in 05 10 15 20 25 30
    do
        ell=5
        for d in 0 1 2
        do
            ./cmm_perf \
                instances/ed${e}/sampled_human_V${V}_ed${e}.fa \
                instances/ed${e}/sampled_human_V${V}_ed${e}.int \
                ${ell} ${d} 0.${e} ${MAX_THREADS} >> ${OUT_RQ2}
        done
    
        ell=8
        for d in 0 1 2 3 4
        do
            ./cmm_perf \
                instances/ed${e}/sampled_human_V${V}_ed${e}.fa \
                instances/ed${e}/sampled_human_V${V}_ed${e}.int \
                ${ell} ${d} 0.${e} ${MAX_THREADS} >> ${OUT_RQ2}
        done
    done
done

## RQ3: memory
#V=100
## binom(8,4)=70, only combination that requests all threads
#ELL=8
#D=4
#ED=05
#for t in 1 2 4 8 16 32 64
#do
#    ./cmm_perf \
#        instances/ed${ED}/sampled_human_V${V}_ed${ED}.fa \
#        instances/ed${ED}/sampled_human_V${V}_ed${ED}.int \
#        ${ELL} ${D} 0.${ED} ${t} >> ${OUT_RQ3}
#
#done
