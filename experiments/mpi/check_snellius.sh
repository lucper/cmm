#!/usr/bin/bash
#
# First job on Snellius: prints what the node is (sockets, NUMA domains, modules), downloads E. coli,
# builds cmm, runs its tests, shows where the threads run, and compares the MPI version with cmm_test.
#
# From the root of the repository:
#   sbatch experiments/mpi/check_snellius.sh
# The results are in experiments/mpi/check_<job id>.log. Other partitions or times: sbatch -p, -t.
#
#SBATCH --job-name=cmm_check
#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-node=8                      # one process per NUMA domain
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00                          # at most 1 hour: the nodes for short jobs
#SBATCH --output=experiments/mpi/check_%j.log    # stdout and stderr

set -u

if [ ! -f src/main.cpp ]
then
    sed -n '3,8p' "$0" >&2
    exit 1
fi

step() { echo; echo "=== $*"; }

# The toolchain the wiki recommends; it lists the others if this one is not there.
if ! module load 2025 foss/2025a
then
    module avail foss
    exit 1
fi

export LC_ALL=C
export OMP_PLACES=cores OMP_PROC_BIND=close
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}  # srun does not take it from #SBATCH

THREADS=${SLURM_CPUS_PER_TASK}
ALL_THREADS=$((SLURM_NTASKS * SLURM_CPUS_PER_TASK))

step "The node"
hostname
lscpu
numactl --hardware
module list
mpirun --version
mpicxx --version

# The options of data/run.sh are all required, but only E. coli is used here; it skips what it has.
step "Dataset"
if ! experiments/data/run.sh -O "511145:Escherichia_coli_K12_MG1655" -s 700 -S 42 -x 511145 -e 05 -v 100
then
    echo "FAILED: no dataset. If this node has no internet, run the data/run.sh command of this script on a login node, then resubmit."
    exit 1
fi
FASTA=experiments/data/real_instances/511145.sequences.s700.clean.fa
EDGES=experiments/data/real_instances/511145.links.s700.clean.int

# Not while experiment jobs run from this clone, as it replaces their programs.
step "Build"
experiments/mpi/build_snellius.sh || exit 1

step "Tests"
make test || echo "FAILED: make test"
make test_mpi || echo "FAILED: make test_mpi"

# Each thread should have its own core, not all those of a process the same one.
step "Where the processes and threads run"
OMP_DISPLAY_AFFINITY=TRUE srun --distribution=block:block --cpu-bind=verbose ./cmm_mpi \
    -s test/instances/readme.fa -i test/instances/readme.int -l 4 -d 2 -f x2 -k 1 -t "${THREADS}"

# The columns of the result lines: experiments/mpi/README.md.
step "One process against ${SLURM_NTASKS} processes, on E. coli"
solutions=$(mktemp -d -p "${TMPDIR:-/tmp}")
for d in 0 1 2 3 4
do
    echo "--- ell = 5, d = ${d}"
    ./experiments/bin/cmm_test "${FASTA}" "${EDGES}" 5 "${d}" 1000 "${ALL_THREADS}" "${solutions}/one" x2 2>/dev/null
    srun --distribution=block:block ./experiments/bin/cmm_test_mpi "${FASTA}" "${EDGES}" 5 "${d}" 1000 "${THREADS}" "${solutions}/mpi" x2 2>/dev/null
    cmp -s "${solutions}/one" "${solutions}/mpi" && echo "same pairs" || echo "FAILED: different pairs"
done
rm -r "${solutions}"
