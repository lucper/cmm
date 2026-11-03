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
On one whole node, it prints the sockets, NUMA domains, and modules, downloads the E. coli dataset with `../data/run.sh`, builds everything, runs `make test` and `make test_mpi`, shows the core of each thread, checks that each of the 8 processes stays in one NUMA domain, and compares `cmm_test_mpi` with `cmm_test` at motif length 5 and $d = 0, \dots, 4$.
After cloning the repository, from its root directory:
```bash
sbatch experiments/mpi/check_snellius.sh
```
If the node has no internet access, the job stops and prints this: run the `data/run.sh` command of the script on a login node, then submit again.
The results are in `experiments/mpi/check_<job id>.log`: any `FAILED` there needs a look, and so does a thread that shares its core.
The module names in the scripts may need adjusting to the machine.

## Layouts of processes and threads

The same program can divide the cores of a node in several ways: for example, on the 128 cores of a `rome` node, 8 processes of 16 threads (one process per NUMA domain, the default), 2 of 64 (one per socket), or 1 of 128 (one per node).
`runme.sh` compares them on whole nodes, end to end, with all the parameters of the experiment at its top, for now `human_ii` at $\ell = 8$ and $d = 5$, on 1, 2, 4, and 8 nodes; $d = 3$ is commented out, for after a first look at the results.
Edit them, and on a login node:
```bash
experiments/mpi/runme.sh
```
It only submits jobs, of which Slurm decides the order, so nothing has to stay running:

1. `build_snellius.sh` builds the programs and unpacks the datasets from SLIDER. It is a job of its own, so that several experiment jobs can run at once from one clone. Set `BUILD=0` while jobs of the same clone are running, as a build replaces their programs; build on the partition the jobs run on, as `-march=native` compiles for the CPU of that node.
2. One job per number of wildcards and of nodes runs `make` in `010-measure/`, for all the layouts and trials in turn. The jobs do not fetch datasets, so run `../data/run.sh` first for the ones from STRING, e.g., on a login node. One more job per number of wildcards runs the current code without MPI, `cmm_test`, on 1 node, to compare with the MPI version and with older code.
3. One job per number of wildcards runs `make` in `020-concat/`, when those have ended, whether they finished or not. It writes the table of results, and logs `FAILED` if the runs disagree on the pairs.

`job.sh` is the job of steps 2 and 3: the Slurm settings and modules, around the given command.
Their logs are in `experiments/mpi/job_<job id>.log`, with the cores of the processes of each run.
The jobs do not stop by themselves if a process spans several NUMA domains, so check a new partition with `check_snellius.sh` first, or read the lines of `grep cpu-bind` in those logs, as soon as a run has started: in a layout with one process per domain, each process must have one block of consecutive cores of its domain, such as `0xffff` and `0xffff0000` for 16 cores.

An older run of `human_ii` with 128 threads, on 1 node, took about 28 minutes at $d = 5$, and 11.3 hours at $d = 3$, with 93 GiB; each job of step 2 does 3 such runs. So the time limits in `runme.sh` are about twice that on 1 node, and leave less margin on 8 nodes.

As in the suites of `../README.md`, `make` decides what to run, one measurement per file; it runs inside the job, as `sbatch` returns before a job has run, and Slurm orders the jobs (`--dependency`).
So after a job fails, or is killed at its time limit, run `runme.sh` again with `BUILD=0`: it skips the runs that have finished, and collects the results again.
Each run is in `010-measure/<dataset>/`, with its solution (`.out`), and all of them are in `020-concat/<dataset>.layouts.l<ell>d<d>.<f>.k<k>.tsv`: each line has the number of nodes, the layout and the trial, then the [result line](#cmm_test_mpi) of `cmm_test_mpi`.
The runs without MPI are there too, as the layout `nompi<threads>`, with the 13 columns of `cmm_test`.
A result file is not redone when the code changes, so move the old results aside first.

## Comparison with ALENEX_AEC

`aec/` is a separate check, to compare the running time and memory of the MPI version with the code of the tag `ALENEX27_AEC`, with `AEC=1` in `runme.sh`; `010-measure/` and `020-concat/` do not depend on it.
`runme.sh` also submits a job that builds the tag in `aec/source/`, one job per number of wildcards that runs its `cmm_test` on 1 node, and one that merges the results, when the collection of step 3 and that run have ended.
This only works with the support function `x2`, the only one of the tag.

The results are in `aec/<dataset>.layouts.l<ell>d<d>.<f>.k<k>.tsv`: the lines of the same file in `020-concat/`, then those of the tag, as the layout `aec<threads>`, with the 13 columns of `cmm_test`.
Its column 10, the prunings, has no match, so it is dropped; it stays in `aec/<dataset>/`, with the solutions of the tag.
These are not compared with the others: the code of the tag breaks ties differently, so check them by hand.
The difference in time also includes what branch 57 changed, the removal of pruning and of the order of the threads in the ties, and not only the refactoring.
