#!/usr/bin/bash
#
# Runs the experiment of the layouts end to end, with all its parameters here. On a login node:
#   experiments/mpi/runme.sh
# It submits the build, one job per number of wildcards and of nodes, and one without MPI per number of
# wildcards, which wait for it, and for each number of wildcards one job that collects the results, which waits
# for those: see README.md. Submitting again skips
# the runs that have finished. Set BUILD=0 while jobs of this clone are running, as a build replaces their
# programs.

set -euo pipefail

BUILD=1
DATASET=human_ii
FASTA=experiments/data/slider_data/human_ii.clean.fa    # from the root of the repository
INTS=experiments/data/slider_data/human_ii.clean.int
ELL=8
SUPPORT=x2
K=10000
LAYOUTS="8x16 2x64 1x128"    # processes per node x threads per process
THREADS=128                  # of the run without MPI
TRIALS=1

# For each number of wildcards, the number of nodes:time limit of its jobs, from the old runs (README.md).
declare -A JOBS=(
    [5]="1:03:00:00 2:02:00:00 4:01:00:00 8:01:00:00"
    [3]="1:60:00:00 2:30:00:00 4:15:00:00 8:08:00:00"
)
# For each of them, the time limit of the run without MPI, on 1 node.
declare -A NO_MPI_TIME=(
    [5]="01:00:00"
    [3]="24:00:00"
)

cd "$(dirname "$0")/../.."

submit() { sbatch --parsable --kill-on-invalid-dep=yes "$@"; }

dependency=()
if [ "${BUILD}" = 1 ]
then
    build=$(submit experiments/mpi/build_snellius.sh)
    echo "build: job ${build}"
    dependency=("--dependency=afterok:${build}")
fi

VARIABLES=(DATASET="${DATASET}" L="${ELL}" F="${SUPPORT}" K="${K}")

for d in "${!JOBS[@]}"
do
    FILES=(FASTA="${PWD}/${FASTA}" INTS="${PWD}/${INTS}")

    id=$(submit "${dependency[@]}" --nodes=1 --exclusive -t "${NO_MPI_TIME[${d}]}" experiments/mpi/job.sh \
             make -k -C experiments/mpi/010-measure no-mpi "${VARIABLES[@]}" D="${d}" "${FILES[@]}" NODES=1 \
             THREADS="${THREADS}" TRIALS="${TRIALS}")
    echo "${DATASET}, ell = ${ELL}, d = ${d}, without MPI, ${THREADS} threads: job ${id}"
    ids=":${id}"

    for job in ${JOBS[${d}]}
    do
        nodes=${job%%:*} time=${job#*:}
        id=$(submit "${dependency[@]}" --nodes="${nodes}" --exclusive -t "${time}" experiments/mpi/job.sh \
                 make -k -C experiments/mpi/010-measure all "${VARIABLES[@]}" D="${d}" "${FILES[@]}" \
                 NODES="${nodes}" LAYOUTS="${LAYOUTS}" TRIALS="${TRIALS}")
        echo "${DATASET}, ell = ${ELL}, d = ${d}, ${nodes} node(s), ${time}: job ${id}"
        ids="${ids}:${id}"
    done
    # afterany: also if a job failed or reached its limit, with the runs that finished.
    id=$(submit --dependency="afterany${ids}" -t 00:10:00 experiments/mpi/job.sh \
             make -C experiments/mpi/020-concat all "${VARIABLES[@]}" D="${d}")
    echo "${DATASET}, ell = ${ELL}, d = ${d}, results: job ${id}"
done
