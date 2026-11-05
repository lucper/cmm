#!/usr/bin/env bash
#
# Fetches the logs of the MPI experiments from Snellius: the check_*, build_* and job_*.log files of the Slurm jobs,
# the timings of the layouts (the .perf and .tsv files of 010-measure/<dataset>/, and the tables of 020-concat/), and
# the results of 025-aec/, without the code of the tag in aec/source/ and the runs still going (*.tmp).
# Files that are gone from Snellius stay here, and a file that changed is fetched again.
#
# From the root of the repository:
#   experiments/mpi/sync-logs.sh
#   HOST=snellius REMOTE=cmm/experiments/mpi experiments/mpi/sync-logs.sh -n    # extra options go to rsync, e.g. a dry run

set -eu
HOST=${HOST:-snellius}
REMOTE=${REMOTE:-cmm/experiments/mpi}
DIR=experiments/mpi

[ -d "$DIR" ] || { echo "Run this from the root of the repository." >&2; exit 1; }

rsync -av --include='*.log' --exclude='*' "$@" "$HOST:$REMOTE/" "$DIR/"
rsync -avm --include='*/' --include='*.perf' --include='*.tsv' --exclude='*' "$@" "$HOST:$REMOTE/010-measure/" "$DIR/010-measure/"
rsync -av --include='*.tsv' --exclude='*' "$@" "$HOST:$REMOTE/020-concat/" "$DIR/020-concat/"
rsync -av --exclude='/source/' --exclude='/Makefile' --exclude='*.tmp' "$@" "$HOST:$REMOTE/025-aec/" "$DIR/025-aec/"
