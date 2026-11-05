#!/usr/bin/bash
#
# Simulates the layouts from the results of experiments/mpi/runme.sh, on your own machine, not in a job: see README.md.
#   experiments/mpi/simulate/runme.sh [name of the results]    # by default human_ii
# The first time, it installs Julia and the packages of the script in .julia/.

set -euo pipefail

DATASET=${1:-human_ii}
L=8 D=5 F=x2 K=10000

cd "$(dirname "$0")"

[ -d .julia ] || make setup
make all DATASET="${DATASET}" L="${L}" D="${D}" F="${F}" K="${K}"
echo "The results are in experiments/mpi/simulate/${DATASET}/."
