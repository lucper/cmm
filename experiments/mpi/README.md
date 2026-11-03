# The MPI version

Here, you will find the job scripts and the experiments of the MPI version, which `cmm_mpi` and `cmm_test_mpi` run on several processes, possibly on several nodes.
The commands are run from the root directory of the repository.

## `cmm_mpi`

`cmm_mpi` runs the search on several processes, possibly on several nodes, each with its own threads.
It needs an MPI implementation with `mpicxx` (e.g., Open MPI), and is built separately, so `cmm` does not need MPI:
```bash
make cmm_mpi
```

It takes the same options as `cmm` and prints the same motif pairs, for any number of processes.
Option `-t` is the number of threads **per process**.
Only the first process reads the input files, so the other nodes do not need them.
For example, on 2 processes with 4 threads each:
```bash
mpirun -np 2 --bind-to none ./cmm_mpi -s example.fa -i example.int -l 5 -d 2 -f x2 -k 5 -t 4
```
`--bind-to none` matters on one machine: with one or two processes, Open MPI may otherwise bind each process to a single core, which all its threads would share.
On a cluster, bind each process to its own cores instead, as the job scheduler allows; for Slurm, see the jobs below.

Each process searches its own share of the wildcard configurations, so there is little to gain from more processes than configurations: a motif length of 5 has at most 55, and a length of 8 has up to 2485.

The tests of the MPI version run with `make test_mpi`.

## `cmm_test_mpi`

The suites of [`../README.md`](../README.md) use `cmm_test`, on one process.
`cmm_test_mpi` is the same program on several processes, built separately so that the suites do not need MPI. From the root directory:
```bash
make -C experiments/src cmm_test_mpi    # builds experiments/bin/cmm_test_mpi; needs mpicxx
```

It takes the arguments of `cmm_test`, with `num_threads` per process, and writes the same solution file.
Its result line has the 13 columns of `cmm_test`, then 5 more:

| Column | Meaning |
|--------|---------|
| 11 | threads spawned, by all the processes |
| 12 | time in ms on the first process, after the merge; it does not include reading the input |
| 13 | peak RAM in KB, the largest of a process |
| 14 | number of processes |
| 15 | peak RAM in KB, the sum over the processes |
| 16, 17, 18 | fastest, median, and slowest search time of a process, in ms |

The peak RAM per thread is column 15 divided by column 11, as column 13 covers only one process.
Columns 16 to 18 show how evenly the configurations are shared: the run waits for the slowest process.

## First run on a cluster

`check_snellius.sh` is a Slurm job that checks a machine before any experiment, here Snellius (partition `rome`).
On one whole node, it prints the sockets, NUMA domains, and modules, downloads the E. coli dataset with `../data/run.sh`, builds everything, runs `make test` and `make test_mpi`, shows the core of each thread, and compares `cmm_test_mpi` with `cmm_test` at motif length 5 and $d = 0, \dots, 4$.
After cloning the repository, from its root directory:
```bash
sbatch experiments/mpi/check_snellius.sh
```
If the node has no internet access, the job stops and prints this: run the `data/run.sh` command of the script on a login node, then submit again.
The results are in `experiments/mpi/check_<job id>.log`: any `FAILED` there needs a look, and so does a thread that shares its core.
The module names at the top of the script may need adjusting to the machine.

## Layouts of processes and threads

The same program can divide the cores of a node in several ways: for example, on the 128 cores of a `rome` node, 8 processes of 16 threads (one process per NUMA domain, the default), 2 of 64 (one per socket), or 1 of 128 (one per node).
`layouts.sh` compares them on whole nodes, as many as `sbatch --nodes` asks for, with one run of each, on E. coli at motif length 8 and the given number of wildcards.
It uses the programs built by `build_snellius.sh`, which is a job of its own, so that several experiment jobs can run at once from one clone: build again only after changing the code, and on the partition they run on, as `-march=native` compiles for the CPU of that node.
From the root directory, for example at $d = 5$, on 1, 2, 4, and 8 nodes:
```bash
job=$(sbatch --parsable experiments/mpi/build_snellius.sh)
for nodes in 1 2 4 8
do
    sbatch --dependency=afterok:${job} --nodes=${nodes} -t 01:00:00 experiments/mpi/layouts.sh 5
done
```
It stops if a process of the default layout does not stay within one NUMA domain, and logs `FAILED` if the layouts disagree on the pairs.

As in the suites of `../README.md`, `make` decides what to run, one measurement per file; it runs inside the job, as `sbatch` returns before a job has run, and Slurm orders the jobs (`--dependency`).
So a job killed at its time limit can be submitted again: it skips the runs that have finished.
Each run is in `010-measure/<dataset>/`, with its solution (`.out`), and all of them are in `020-concat/<dataset>.layouts.l<ell>d<d>.k<k>.tsv`: each line has the number of nodes, the layout and the trial, then the [result line](#cmm_test_mpi) of `cmm_test_mpi`.
A result file is not redone when the code changes, so move the old results aside first.
