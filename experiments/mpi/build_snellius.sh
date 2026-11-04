#!/usr/bin/bash
#
# Builds cmm, cmm_mpi and the programs of experiments/ on a Snellius node, and unpacks the datasets from
# SLIDER, for the experiment jobs, which only use them, so that several can run at once. Build again after changing the code, before any of
# them runs, and on the partition they run on, as -march=native compiles for the CPU of this node.
#
# Submitted by runme.sh, or from the root of the repository:
#   sbatch experiments/mpi/build_snellius.sh
#
#SBATCH --job-name=cmm_build
#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:30:00
#SBATCH --output=experiments/mpi/build_%j.log

module load 2025 foss/2025a || exit 1

[ -d experiments/data/slider_data ] || tar -xzf experiments/data/slider_data.tar.gz -C experiments/data || exit 1

git log -1 --format='Building commit %h'
make clean && make -C experiments/src clean
make -j "$(nproc)" cmm cmm_mpi && make -j "$(nproc)" -C experiments/src all cmm_test_mpi
