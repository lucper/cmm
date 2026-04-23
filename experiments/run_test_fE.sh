#!/usr/bin/bash

NODE_FILE=$1
EDGE_FILE=$2
L=$3
D=$4
K=$5

MAIN_PROG_CMD="../build/cmm -n $NODE_FILE -e $EDGE_FILE -l $L -d $D -k $K -f E"
TEST_PROG_CMD="./main_algo.py $L $D $K $NODE_FILE $EDGE_FILE"

if diff <($MAIN_PROG_CMD | ./format_output.sh) <($TEST_PROG_CMD | ./format_output.sh) > /dev/null
then
    echo "TEST PASSED."
else
    echo "TEST FAILED"
fi
