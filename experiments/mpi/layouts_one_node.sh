#!/usr/bin/bash
#
# Compares layouts of processes and threads on one whole node: cmm_test_mpi on E. coli, at ell = 8, in
# turn. Each run is one file of 010-measure/, skipped if it exists, so a job killed at its time limit
# can be submitted again. All the runs are then in 020-concat/<dataset>.layouts.l<ell>d<d>.k<k>.tsv,
# each line with the layout and the trial; its columns are in experiments/README.md.
#
# From the root of the repository; expect about 3 hours at d = 3, and 20 minutes at d = 5, which can ask
# for -t 01:00:00 to run on the nodes for short jobs:
#   sbatch experiments/mpi/layouts_one_node.sh <d>
# It uses the programs of build_snellius.sh, built on the same partition. Other times: sbatch -t.
#
#SBATCH --job-name=cmm_layouts
#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-node=8                         # the default layout: one process per NUMA domain
#SBATCH --cpus-per-task=16
#SBATCH --time=08:00:00
#SBATCH --output=experiments/mpi/layouts_%j.log     # stdout and stderr

set -u

if [ $# -ne 1 ] || [ ! -x experiments/bin/cmm_test_mpi ]
then
    sed -n '3,11p' "$0" >&2
    exit 1
fi
D=$1
DATASET=511145 ELL=8 K=10000 TRIALS=1
LAYOUTS="8x16 2x64 1x128"    # processes x threads per process

module load 2025 foss/2025a || exit 1
export LC_ALL=C OMP_PLACES=cores OMP_PROC_BIND=close

experiments/data/run.sh -O "511145:Escherichia_coli_K12_MG1655" -s 700 -S 42 -x 511145 -e 05 -v 100 || exit 1
FASTA=$(realpath "experiments/data/real_instances/${DATASET}.sequences.s700.clean.fa")
INTS=$(realpath "experiments/data/real_instances/${DATASET}.links.s700.clean.int")

# Where the processes of a layout run.
where() {
    srun --nodes=1 --ntasks="${1%x*}" --ntasks-per-node="${1%x*}" --cpus-per-task="${1#*x}" \
         --distribution=block:block numactl --show
}
for layout in ${LAYOUTS}
do
    echo "--- where ${layout} runs"
    where "${layout}" | grep -E '^(nodebind|physcpubind)' | sort
done
# Hours of runs would not measure the default layout if a process spanned several NUMA domains.
if [ "$(where 8x16 | awk '/^nodebind/ { print NF - 1 }' | sort -u)" != 1 ]
then
    echo "FAILED: in 8x16, a process is not in exactly one NUMA domain"
    exit 1
fi

NAME=${DATASET}.layouts.l${ELL}d${D}.k${K}
inputs="" first=""
for trial in $(seq 1 "${TRIALS}")
do
    for layout in ${LAYOUTS}
    do
        echo "--- ${layout}, trial ${trial}"
        make -s -C experiments/mpi/010-measure all DATASET="${DATASET}" FASTA="${FASTA}" INTS="${INTS}" \
             L="${ELL}" D="${D}" K="${K}" LAYOUT="${layout}" TRIAL="${trial}" || { echo "FAILED: the run"; continue; }

        run=${DATASET}/${NAME}.${layout}.trial${trial}
        inputs="${inputs} ../010-measure/${run}.perf"
        out=experiments/mpi/010-measure/${run}.out
        first=${first:-${out}}
        cmp -s "${first}" "${out}" || echo "FAILED: different pairs from ${first}"
    done
done

if [ -n "${inputs}" ]
then
    make -s -C experiments/mpi/020-concat all OUT="${NAME}.tsv" INPUTS="${inputs}"
    echo "Done: experiments/mpi/020-concat/${NAME}.tsv"
fi
