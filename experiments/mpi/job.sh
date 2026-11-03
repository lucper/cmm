#!/usr/bin/bash
#
# A Snellius job that runs the given command with the modules of the build: runme.sh submits it, and states
# what depends on the run, such as the nodes and the time, as options of sbatch.
#
#SBATCH --job-name=cmm
#SBATCH --partition=rome
#SBATCH --output=experiments/mpi/job_%j.log    # stdout and stderr

module load 2025 foss/2025a || exit 1
export LC_ALL=C OMP_PLACES=cores OMP_PROC_BIND=close

"$@"
