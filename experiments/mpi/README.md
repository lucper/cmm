# MPI Extension

The experiment runs on a cluster with Slurm, here Snellius (partition `rome`).
All the commands are run from the root directory of the repository.

## Prerequisites

1. An MPI implementation with `mpicxx` (e.g., Open MPI).
   The scripts load the modules `2025` and `foss/2025a` of Snellius (GCC 14.2.0 and Open MPI 5.0.7); `job.sh`, `build_snellius.sh`, and `check_snellius.sh` each need changing on another machine, or when these are retired.
2. The datasets: `build_snellius.sh` unpacks the ones from SLIDER, and the ones from STRING must be fetched first with `experiments/data/run.sh` (see `experiments/README.md`), e.g., on a login node.
3. On a new partition, a first run of `check_snellius.sh` (see [Checking a machine](#checking-a-machine)).

## Usage

On a login node:
```bash
experiments/mpi/runme.sh
```

It shows the parameters of a preset and the jobs they make, lets you change them, and submits the jobs.
It only submits jobs, of which Slurm decides the order, so nothing has to stay running.
The keys: `p` the next preset, `1` to `9` and `a` to `d` change a parameter, ⏎ submits (after asking), `n` prints the `sbatch` commands instead, and `q` quits.
Next to each job, the screen says how many of its runs have finished, which the job skips, and it warns when a build would replace the programs of jobs of the same clone that are queued or running.

| Preset | Runs |
|--------|------|
| `main` (the default) | the layouts on 1, 2, 4, and 8 nodes, the run without MPI, and the code of the tag `ALENEX27_AEC` |
| `timers` | the layouts on 1 node and the run without MPI, under the name `human_ii_timers`, for the time of each configuration (see `simulate_tasks.jl`) |

Without the screen:
```bash
experiments/mpi/runme.sh --yes                       # submit the preset main as it is
experiments/mpi/runme.sh --preset timers --yes       # or another preset
experiments/mpi/runme.sh --yes --dry-run             # print the sbatch commands, and submit nothing
```

The parameters, set by the presets at the top of `runme.sh`:

| Variable | Meaning |
|----------|---------|
| `BUILD` | `1` to submit the build; turn it off while jobs of the same clone are running, as a build replaces their programs |
| `DATASET` | name of the results (used in every output filename); change it to run again the runs that have finished |
| `FASTA` | FASTA file with protein sequences |
| `INTS` | interactions file (`u v` per line) |
| `ELL` | motif length $\ell$ |
| `SUPPORT` | support function, `E` or `x2` |
| `K` | number $K$ of top motif pairs |
| `LAYOUTS` | layouts of the MPI version, as `<processes per node>x<threads per process>` (e.g., `8x16`) |
| `THREADS` | threads of the run without MPI |
| `TRIALS` | number of trials per configuration |
| `AEC` | `1` to also run the code of the tag `ALENEX27_AEC`, which only has the support function `x2` |
| `AEC_THREADS` | threads of the run of the tag |
| `WILDCARDS` | the numbers $d$ of wildcards to run |
| `JOBS` | for each number $d$ of wildcards, the jobs of the MPI version, as `<nodes>:<time limit>` |
| `NO_MPI_TIME` | for each number $d$ of wildcards, the time limit of the runs on 1 node, without MPI |

`runme.sh` submits these jobs:

| Job | Runs | How many | Resources | Starts |
|-----|------|----------|-----------|--------|
| Build | `build_snellius.sh` | 1 | 16 cores | at once, if `BUILD=1` |
| MPI runs | `make -k -C experiments/mpi/010-measure all` | one per $d$ and number of nodes | whole nodes | if the build succeeded |
| Run without MPI | `make -k -C experiments/mpi/010-measure no-mpi` | one per $d$ | 1 whole node | if the build succeeded |
| Table of results | `make -C experiments/mpi/020-concat all` | one per $d$ | smallest allocation | when the two above have ended, whether they finished or not |
| Build of the tag | `make -C experiments/mpi/aec build` | 1, if `AEC=1` | 16 cores | at once |
| Run of the tag | `make -k -C experiments/mpi/aec run` | one per $d$, if `AEC=1` | 1 whole node | if its build succeeded |
| Merged table | `make -C experiments/mpi/aec merge` | one per $d$, if `AEC=1` | smallest allocation | when the table of results and the run of the tag have ended |

All of them but the build run through `job.sh`, which sets the Slurm options and loads the modules.

The experiment is a sequence of numbered `make` stages, each consuming the previous stage's output, and `aec/` is a separate check that reads the result of the last one:

| Directory | Produces |
|-----------|----------|
| `010-measure/` | One measurement per file, for each layout and trial, and for the run without MPI |
| `020-concat/` | One TSV per dataset, $\ell$, $d$, support function, and $K$, with the measurements of all the jobs, and a check that they found the same pairs |
| `aec/` | The same TSV, with the measurements of the code of the tag |

`make` runs inside the jobs, as `sbatch` returns before a job has run, and Slurm orders the jobs (`--dependency`).
After a job fails or is killed at its time limit, run `runme.sh` again with the build off: it skips the runs that have finished, and collects the results again.
A measurement is not redone when the code changes, so move the old results aside first.

### Layouts

The same program can divide the cores of a node in several ways, which differ only in how it is launched.
On the 128 cores of a `rome` node:

| Layout | Processes per node | Threads per process | One process per |
|--------|--------------------|---------------------|-----------------|
| `8x16` | 8 | 16 | NUMA domain |
| `2x64` | 2 | 64 | socket |
| `1x128` | 1 | 128 | node |

Slurm assigns consecutive cores to each process (`--distribution=block:block`), which in these layouts places each process within its NUMA domain or socket.
The jobs do not check this.
So, on a new partition, run `check_snellius.sh` first, which logs `FAILED` if a process spans several NUMA domains.
You can also check a running experiment with `grep cpu-bind experiments/mpi/job_*.log`: with one process per domain of 16 cores, each process should show one block of 16 consecutive cores, such as `0xffff` or `0xffff0000`.

### Checking a machine

`check_snellius.sh` is a job that checks a machine before any experiment:
```bash
sbatch experiments/mpi/check_snellius.sh
```

On one whole node, in this order:

| Step | Shows |
|------|-------|
| The node | sockets, NUMA domains, modules, and the versions of MPI and the compiler |
| Dataset | the download of the E. coli dataset with `experiments/data/run.sh` |
| Build | the build of everything |
| Tests | `make test` and `make test_mpi` |
| Where the processes and threads run | the core of each thread |
| One process in each NUMA domain | whether each of the 8 processes stays in one domain |
| One process against 8 processes | `cmm_test_mpi` against `cmm_test`, at motif length 5 and $d = 0, \dots, 4$ |

The results are in `experiments/mpi/check_<job id>.log`.
Any `FAILED` there needs a look, and so does a thread that shares its core.

## Files

| File / dir | Role |
|------------|------|
| `runme.sh` | Driver: holds the presets of the parameters, shows and edits them, and submits all the jobs. |
| `job.sh` | A job: the Slurm settings and modules, around the given command. |
| `build_snellius.sh` | The job that builds the programs, and unpacks the datasets from SLIDER. |
| `check_snellius.sh` | The job that checks a machine. |
| `<NNN>-*/` | The numbered pipeline stages, each with its own `Makefile`. |
| `aec/` | The check of the code of the tag, with its own `Makefile`. |
| `aec/source/` | The code of the tag, built by `make build` (generated); remove it to build again on another partition. |
| `010-measure/<dataset>/`, `aec/<dataset>/` | The measurements and solutions (generated). |
| `job_<job id>.log`, `build_<job id>.log`, `check_<job id>.log` | The logs of the jobs (generated), with the cores of the processes of each run. |

## Output format

Filenames encode the full parameter set, so outputs from different parameters never collide.
A representative name, of a measurement in `010-measure/`, with its solution in a `.out` file of the same name:

```
<dataset>/<dataset>.layouts.l<L>d<D>.<F>.k<K>.n<nodes>.<layout>.trial<i>.perf
```

where `<layout>` is a layout, `nompi<threads>` for the run without MPI, or `aec<threads>` for the run of the tag.
The TSVs are `020-concat/<dataset>.layouts.l<L>d<D>.<F>.k<K>.tsv` and, with the same name, in `aec/`.

Each row of the TSVs is one measurement (the one-line output of `cmm_test_mpi` or `cmm_test`) with the number of nodes, the layout, and the trial prepended, tab-separated:

| # | Column | | # | Column | | # | Column |
|---|--------|-|---|--------|-|---|--------|
| 1 | `nodes` | | 8 | `ell` | | 15 | `time_ms` |
| 2 | `layout` | | 9 | `d` | | 16 | `peak_ram_kb` |
| 3 | `trial` | | 10 | `max_rank` | | 17 | `num_processes` |
| 4 | `instance` | | 11 | `min_rank` | | 18 | `sum_peak_ram_kb` |
| 5 | `edge_density` | | 12 | `avg_rank` | | 19 | `min_search_ms` |
| 6 | `V` | | 13 | `num_threads_requested` | | 20 | `median_search_ms` |
| 7 | `N` | | 14 | `num_threads_spawned` | | 21 | `max_search_ms` |

Columns 17 to 21 are only in the rows of `cmm_test_mpi`, and in these rows:
- `num_threads_requested` (13) is per process, and `num_threads_spawned` (14) is over all the processes;
- `time_ms` (15) is on the first process, after the merge, and does not include reading the input;
- `peak_ram_kb` (16) is the largest of a process, and `sum_peak_ram_kb` (18) the sum over the processes, so the peak RAM per thread is (18) divided by (14);
- `min_search_ms`, `median_search_ms`, and `max_search_ms` (19 to 21) are the search times of the processes, which show how evenly the configurations are shared, as the run waits for the slowest process; with an even number of processes, the median is the upper one.

In the rows of the tag, its column of prunings has no match, so it is dropped; it stays in its `.perf` file.

## Programs

`cmm_mpi` runs the search on several processes, possibly on several nodes, each with its own threads.
It is built separately, so `cmm` does not need MPI:
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

Each process searches its own share of the wildcard configurations, so there is little to gain from more processes than configurations: a motif length of 5 has at most 55, and a length of 8 has up to 2485.

`cmm_test_mpi` is `cmm_test` on several processes, built separately so that the suites of `experiments/README.md` do not need MPI:
```bash
make -C experiments/src cmm_test_mpi    # builds experiments/bin/cmm_test_mpi
```
It takes the arguments of `cmm_test`, with `num_threads` per process, and writes the same solution file.

The tests of the MPI version run with `make test_mpi`.
