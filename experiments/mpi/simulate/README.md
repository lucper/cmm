# The simulation of the layouts

This directory simulates how the search would run on more nodes, for testing and benchmarking.
It is not part of the main pipeline (`experiments/mpi/runme.sh`), and no job runs it: it runs on your own machine, on the results of that pipeline.

## Layouts

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

To measure the other layouts, add them on the screen of `experiments/mpi/runme.sh` (key `a`).

## Simulating the layouts

`simulate_tasks.jl` simulates the search on more nodes from the time of each configuration on 1 node, to show what limits the speedup.
It reads the results of `010-measure/` and `020-concat/` of `experiments/mpi/runme.sh`, of each layout measured:
```bash
experiments/mpi/simulate/runme.sh                     # the results named human_ii, L=8 D=5 F=x2 K=10000
experiments/mpi/simulate/runme.sh human_ii_round3     # results of another name
```

It installs Julia the first time (`make setup`), and then runs `make all` of the `Makefile` here, whose variables you can also give yourself:
```bash
make -C experiments/mpi/simulate all DATASET=human_ii L=8 D=5 F=x2 K=10000
```

`make setup` installs Julia and the packages of `simulate_tasks.jl` in `simulate/.julia/`, about 2 GB on Linux x86-64, and nothing elsewhere; to remove them, delete that directory.
On another system, give another Julia with `JULIA=<program>`, and install the packages of `Project.toml` yourself.

If the results are still on Snellius, the helper `experiments/mpi/sync-logs.sh` copies them to the same directories on your machine (see `../README.md`).

The results are in `simulate/<dataset>/`, named as the TSV of `020-concat/`:

| File | Holds |
|------|-------|
| `<name>.txt` | The printed results of all the layouts |
| `<name>-<layout>-speedup.pdf`, `.png` | The measured against the simulated speedup |
| `<name>-<layout>.pdf`, `.png` | The same, above the times of the configurations, the simulated time, and its spread over the trials |

The simulation assumes that:
- each layout has its times on 1 node, or it is skipped without an error, and the table has the numbers of nodes to simulate, so a dataset with only 1 node simulates only that;
- the times come from one trial (`TRIAL`, by default 1), while the measured values of the table are the mean of all its trials, so they compare best with one trial per run;
- a configuration takes the same time on any number of nodes, and every node has the speeds of the cores and warm-ups of the measured node;
- the code still splits the configurations into blocks and schedules them as the script does, by hand, in `simulate_tasks.jl`, and the columns of the `.perf` files and the TSV stay in place.

The script runs on a fixed number of threads, so with the same version of Julia it gives the same results.
