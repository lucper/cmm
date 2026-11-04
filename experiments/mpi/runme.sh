#!/usr/bin/bash
#
# Runs the experiment of the layouts end to end. On a login node:
#   experiments/mpi/runme.sh                    # choose and edit what to run, then submit
#   experiments/mpi/runme.sh --yes              # submit the preset main as it is, without asking
#   experiments/mpi/runme.sh --preset timers    # start from another preset
#   experiments/mpi/runme.sh --dry-run          # print the sbatch commands, and submit nothing
# It submits the build, one job per number of wildcards and of nodes, and one without MPI per number of
# wildcards, which wait for it, and for each number of wildcards one job that collects the results, which waits
# for those: see README.md. Submitting again skips the runs that have finished. Turn the build off while jobs
# of this clone are running, as a build replaces their programs.
#
# The presets, of which the screen shows the parameters:
# - main: the layouts on 1, 2, 4 and 8 nodes, with the code of the tag ALENEX27_AEC (aec/);
# - timers: the layouts on 1 node only, under their own name, for the time of each configuration that the code
#   prints with the timers (see 010-measure/Makefile and simulate_tasks.jl).

set -euo pipefail

# The parameters, which a preset sets, and the screen edits.
declare -A JOBS NO_MPI_TIME
PRESETS=(main timers)

# The layouts on 1, 2, 4 and 8 nodes, with ALENEX27_AEC.
preset_main() {
    PRESET=main
    DESCRIPTION="the layouts on 1, 2, 4 and 8 nodes, with ALENEX27_AEC"
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
    AEC=1                        # also run the code of the tag ALENEX27_AEC, and merge it with the results (aec/)
    AEC_THREADS=128
    WILDCARDS="5"                # the numbers of wildcards to run, of those of JOBS
    # For each number of wildcards, the number of nodes:time limit of its jobs, from the old runs (README.md).
    JOBS=([5]="1:03:00:00 2:02:00:00 4:01:00:00 8:01:00:00"
          [3]="1:60:00:00 2:30:00:00 4:15:00:00 8:08:00:00")
    # For each of them, the time limit of the run without MPI, on 1 node.
    NO_MPI_TIME=([5]="01:00:00" [3]="24:00:00")
    EDITED=0
}

# The layouts on 1 node, under their own name, to time each configuration.
preset_timers() {
    preset_main
    PRESET=timers
    DESCRIPTION="the layouts on 1 node, to time each configuration"
    DATASET=human_ii_timers
    AEC=0
    JOBS=([5]="1:03:00:00" [3]="1:60:00:00")
}

usage() {
    sed -n '3,7s/^# \{0,1\}//p' "$0"
    echo "Presets: ${PRESETS[*]} (default: ${PRESETS[0]})."
}

DRY_RUN=0
ASK=1
preset=${PRESETS[0]}
while [ $# -gt 0 ]
do
    case "$1" in
        -y|--yes) ASK=0 ;;
        -n|--dry-run) DRY_RUN=1 ;;
        -p|--preset) preset=${2:?--preset needs a name}; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
    shift
done
[[ " ${PRESETS[*]} " == *" ${preset} "* ]] || { echo "Unknown preset: ${preset}" >&2; usage >&2; exit 2; }
"preset_${preset}"

cd "$(dirname "$0")/../.."

# ------------------------------------------------------------------------------------------------------------------
# What will run: the jobs, and how many of their runs have finished, whose results are there already.

# The number of finished runs among the .perf files given.
count_done() {
    local n=0 f
    for f in "$@"; do [ -e "$f" ] && n=$((n + 1)); done
    echo "${n}"
}

# The .perf files of the runs of the layouts of d on nodes nodes, and of the run without MPI.
layout_perfs() {
    local d=$1 nodes=$2 trial layout
    for trial in $(seq 1 "${TRIALS}"); do
        for layout in ${LAYOUTS}; do
            echo "experiments/mpi/010-measure/${DATASET}/${DATASET}.layouts.l${ELL}d${d}.${SUPPORT}.k${K}.n${nodes}.${layout}.trial${trial}.perf"
        done
    done
}
no_mpi_perfs() {
    local d=$1 trial
    for trial in $(seq 1 "${TRIALS}"); do
        echo "experiments/mpi/010-measure/${DATASET}/${DATASET}.layouts.l${ELL}d${d}.${SUPPORT}.k${K}.n1.nompi${THREADS}.trial${trial}.perf"
    done
}
aec_perfs() {
    local d=$1 trial
    for trial in $(seq 1 "${TRIALS}"); do
        echo "experiments/mpi/aec/${DATASET}/${DATASET}.layouts.l${ELL}d${d}.${SUPPORT}.k${K}.n1.aec${AEC_THREADS}.trial${trial}.perf"
    done
}

# "k of n runs done", for the runs of a job, if any have finished: the job skips them.
done_note() {
    local -a perfs
    mapfile -t perfs < <("$@")
    local done
    done=$(count_done "${perfs[@]}")
    if [ "${done}" = "${#perfs[@]}" ]; then echo "${GREEN}all $(plural "${done}" run) done${RESET}"
    elif [ "${done}" != 0 ]; then echo "${done} of $(plural "${#perfs[@]}" run) done"
    fi
}

plural() { [ "$1" = 1 ] && echo "$1 $2" || echo "$1 $2s"; }

# Prints the jobs that would be submitted, one per line, and sets NUM_JOBS.
plan() {
    NUM_JOBS=0
    local d job nodes time note
    if [ "${BUILD}" = 1 ]; then
        printf '  %-40s %s\n' "build" "build_snellius.sh"
        NUM_JOBS=$((NUM_JOBS + 1))
    fi
    if [ "${AEC}" = 1 ]; then
        printf '  %-40s %s\n' "build of ALENEX27_AEC" "00:30:00"
        NUM_JOBS=$((NUM_JOBS + 1))
    fi
    for d in ${WILDCARDS}; do
        note=$(done_note no_mpi_perfs "${d}")
        printf '  %-40s %-12s %s\n' "d = ${d}: without MPI, ${THREADS} threads" "${NO_MPI_TIME[${d}]}" "${note}"
        NUM_JOBS=$((NUM_JOBS + 1))
        for job in ${JOBS[${d}]}; do
            nodes=${job%%:*} time=${job#*:}
            note=$(done_note layout_perfs "${d}" "${nodes}")
            printf '  %-40s %-12s %s\n' "       $(plural "${nodes}" node), ${LAYOUTS// /, }" "${time}" "${note}"
            NUM_JOBS=$((NUM_JOBS + 1))
        done
        printf '  %-40s %s\n' "       collect the results (020-concat)" "00:10:00"
        NUM_JOBS=$((NUM_JOBS + 1))
        if [ "${AEC}" = 1 ]; then
            note=$(done_note aec_perfs "${d}")
            printf '  %-40s %-12s %s\n' "       ALENEX27_AEC, ${AEC_THREADS} threads" "${NO_MPI_TIME[${d}]}" "${note}"
            printf '  %-40s %s\n' "       merge with ALENEX27_AEC (aec/)" "00:10:00"
            NUM_JOBS=$((NUM_JOBS + 2))
        fi
    done
}

# ------------------------------------------------------------------------------------------------------------------
# The screen: the parameters, the jobs they make, and the keys to change them.

if [ -t 1 ]; then
    BOLD=$'\e[1m' DIM=$'\e[2m' RED=$'\e[31m' YELLOW=$'\e[33m' GREEN=$'\e[32m' RESET=$'\e[0m'
else
    BOLD="" DIM="" RED="" YELLOW="" GREEN="" RESET=""
fi

yes_no() { [ "$1" = 1 ] && echo "yes" || echo "no"; }
exists_note() { [ -e "$1" ] || echo " ${RED}(missing)${RESET}"; }

# The jobs of this clone queued or running, as Slurm reports them: a build would replace their programs.
running_jobs() {
    command -v squeue > /dev/null || { echo 0; return; }
    squeue --me --noheader --format=%Z 2> /dev/null | grep -c -x -F "${PWD}" || true
}

field() { printf '  %s%2s%s  %-22s %s\n' "${BOLD}" "$1" "${RESET}" "$2" "$3"; }

draw() {
    printf '\e[H\e[2J'
    local edited="" d jobs=""
    [ "${EDITED}" = 1 ] && edited=" ${YELLOW}(edited)${RESET}"
    for d in ${WILDCARDS}; do jobs+="d=${d}: ${JOBS[${d}]}; "; done
    echo "${BOLD}The experiment of the layouts${RESET}  ${DIM}experiments/mpi/runme.sh${RESET}"
    echo
    field p "preset" "${BOLD}${PRESET}${RESET}${edited}  ${DIM}${DESCRIPTION}${RESET}"
    echo
    field 1 "build first" "$(yes_no "${BUILD}")"
    field 2 "name of the results" "${DATASET}"
    field 3 "sequences" "${FASTA}$(exists_note "${FASTA}")"
    field 4 "interactions" "${INTS}$(exists_note "${INTS}")"
    field 5 "motif length" "${ELL}"
    field 6 "wildcards" "${WILDCARDS}"
    field 7 "nodes:time limit" "${jobs%; }"
    field 8 "support function" "${SUPPORT}"
    field 9 "k" "${K}"
    field a "layouts" "${LAYOUTS}"
    field b "trials" "${TRIALS}"
    field c "threads without MPI" "${THREADS}"
    field d "ALENEX27_AEC" "$(yes_no "${AEC}")$([ "${AEC}" = 1 ] && echo ", ${AEC_THREADS} threads")"
    echo
    local lines
    lines=$(plan)
    plan > /dev/null
    echo "${BOLD}$(plural "${NUM_JOBS}" job)${RESET}, of which Slurm decides the order:"
    echo "${lines}"
    if [[ "${lines}" == *" done"* ]]; then
        echo "${DIM}The jobs skip the runs done; to run them again, change the name of the results (2).${RESET}"
    fi
    local running
    running=$(running_jobs)
    if [ "${BUILD}" = 1 ] && [ "${running}" != 0 ]; then
        echo
        echo "${YELLOW}$(plural "${running}" job) of this clone queued or running: the build would replace their programs.${RESET}"
    fi
    [ -n "${MESSAGE:-}" ] && { echo; echo "${MESSAGE}"; MESSAGE=""; }
    echo
    echo "${DIM}p${RESET} preset  ${DIM}1-9 a-d${RESET} change  ${DIM}⏎${RESET} submit  ${DIM}n${RESET} dry run  ${DIM}q${RESET} quit"
}

# Reads a new value for a parameter, starting from the current one; keeps it if check fails.
edit() {
    local name=$1 label=$2 pattern=$3 value
    printf '\e[?25h'
    read -r -e -i "${!name}" -p "  ${label}: " value || return 0
    value=$(echo "${value}" | xargs)  # trim, and squeeze the spaces
    if [[ "${value}" =~ ${pattern} ]]; then
        [ "${value}" != "${!name}" ] && EDITED=1
        printf -v "${name}" '%s' "${value}"
    else
        MESSAGE="${RED}Not valid for ${label}: '${value}'${RESET}"
    fi
}

# Reads the jobs and the time limit without MPI of each number of wildcards.
edit_jobs() {
    local d value
    for d in ${WILDCARDS}; do
        value=${JOBS[${d}]:-1:01:00:00}
        read -r -e -i "${value}" -p "  d = ${d}, nodes:time limit of each job: " value || return 0
        value=$(echo "${value}" | xargs)
        if [[ "${value}" =~ ^([0-9]+:[0-9:-]+)( [0-9]+:[0-9:-]+)*$ ]]; then
            [ "${value}" != "${JOBS[${d}]:-}" ] && EDITED=1
            JOBS[${d}]=${value}
        else
            MESSAGE="${RED}Not a list of nodes:time limit: '${value}'${RESET}"
            return 0
        fi
        value=${NO_MPI_TIME[${d}]:-01:00:00}
        read -r -e -i "${value}" -p "  d = ${d}, time limit without MPI: " value || return 0
        value=$(echo "${value}" | xargs)
        if [[ "${value}" =~ ^[0-9:-]+$ ]]; then
            [ "${value}" != "${NO_MPI_TIME[${d}]:-}" ] && EDITED=1
            NO_MPI_TIME[${d}]=${value}
        else
            MESSAGE="${RED}Not a time limit: '${value}'${RESET}"
        fi
    done
}

# Every number of wildcards needs its jobs: asks for those of the new ones.
check_wildcards() {
    local d
    for d in ${WILDCARDS}; do
        [ -n "${JOBS[${d}]:-}" ] && [ -n "${NO_MPI_TIME[${d}]:-}" ] || { edit_jobs; return; }
    done
}

next_preset() {
    local i
    for i in "${!PRESETS[@]}"; do
        if [ "${PRESETS[$i]}" = "${PRESET}" ]; then
            "preset_${PRESETS[$(( (i + 1) % ${#PRESETS[@]} ))]}"
            return
        fi
    done
}

toggle() { printf -v "$1" '%s' $(( 1 - ${!1} )); EDITED=1; }

INT='^[1-9][0-9]*$'
screen() {
    trap 'printf "\e[?25h"' EXIT
    local key rest
    while true; do
        draw
        printf '\e[?25l'
        IFS= read -r -s -n 1 key || exit 1
        # An arrow or another key with an escape sequence: drop the rest of it.
        [ "${key}" = $'\e' ] && read -r -s -t 0.05 -n 5 rest || true
        printf '\e[?25h'
        case "${key}" in
            p) next_preset ;;
            1) toggle BUILD ;;
            2) edit DATASET "name of the results" '^[A-Za-z0-9_.-]+$' ;;
            3) edit FASTA "sequences" '^[^ ]+$' ;;
            4) edit INTS "interactions" '^[^ ]+$' ;;
            5) edit ELL "motif length" "${INT}" ;;
            6) edit WILDCARDS "wildcards" '^[0-9]+( [0-9]+)*$'; check_wildcards ;;
            7) edit_jobs ;;
            8) edit SUPPORT "support function" '^[A-Za-z0-9_]+$' ;;
            9) edit K "k" "${INT}" ;;
            a) edit LAYOUTS "layouts" '^[0-9]+x[0-9]+( [0-9]+x[0-9]+)*$' ;;
            b) edit TRIALS "trials" "${INT}" ;;
            c) edit THREADS "threads without MPI" "${INT}" ;;
            d) toggle AEC; [ "${AEC}" = 1 ] && edit AEC_THREADS "threads of ALENEX27_AEC" "${INT}" ;;
            n) DRY_RUN=1; return ;;
            q) echo; exit 0 ;;
            "")
                plan > /dev/null
                read -r -p "  Submit $(plural "${NUM_JOBS}" job) to Slurm? [y/N] " key || exit 1
                [[ "${key}" == [yY]* ]] && return
                ;;
        esac
    done
}

if [ "${ASK}" = 1 ]; then
    [ -t 0 ] && [ -t 1 ] || { echo "Not a terminal: pass --yes to submit the preset ${PRESET} as it is." >&2; exit 2; }
    screen
    printf '\e[H\e[2J'
fi

# ------------------------------------------------------------------------------------------------------------------
# The jobs.

for d in ${WILDCARDS}; do
    [ -n "${JOBS[${d}]:-}" ] && [ -n "${NO_MPI_TIME[${d}]:-}" ] || { echo "No jobs for d = ${d}" >&2; exit 2; }
done

echo "${BOLD}Preset ${PRESET}$([ "${EDITED}" = 1 ] && echo ", edited")${RESET}: ${DATASET}, ell = ${ELL}, d = ${WILDCARDS}, layouts ${LAYOUTS}"
[ "${DRY_RUN}" = 1 ] && echo "${YELLOW}Dry run: the sbatch commands, with made-up job ids; nothing is submitted.${RESET}"

# The made-up job ids of a dry run count in a file, as submit runs in a subshell, $(submit ...).
if [ "${DRY_RUN}" = 1 ]; then
    fake_ids=$(mktemp)
    trap 'rm -f "${fake_ids}"' EXIT
fi
submit() {
    if [ "${DRY_RUN}" = 1 ]; then
        echo >> "${fake_ids}"
        printf '%s\n' "${DIM}sbatch --parsable --kill-on-invalid-dep=yes $(printf '%q ' "$@")${RESET}" >&2
        echo $((1000 + $(wc -l < "${fake_ids}")))
    else
        sbatch --parsable --kill-on-invalid-dep=yes "$@"
    fi
}

dependency=()
if [ "${BUILD}" = 1 ]
then
    build=$(submit experiments/mpi/build_snellius.sh)
    echo "build: job ${build}"
    dependency=("--dependency=afterok:${build}")
fi

VARIABLES=(DATASET="${DATASET}" L="${ELL}" F="${SUPPORT}" K="${K}")

aec_build=""
if [ "${AEC}" = 1 ]
then
    aec_build=$(submit --ntasks=1 --cpus-per-task=16 -t 00:30:00 experiments/mpi/job.sh \
                    make -C experiments/mpi/aec build)
    echo "build of ALENEX27_AEC: job ${aec_build}"
fi

for d in ${WILDCARDS}
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
    collect=$(submit --dependency="afterany${ids}" -t 00:10:00 experiments/mpi/job.sh \
                  make -C experiments/mpi/020-concat all "${VARIABLES[@]}" D="${d}")
    echo "${DATASET}, ell = ${ELL}, d = ${d}, results: job ${collect}"

    # The code of the tag ALENEX27_AEC, apart from the above; the merge waits for both.
    if [ "${AEC}" = 1 ]
    then
        aec=$(submit --dependency="afterok:${aec_build}" --nodes=1 --exclusive -t "${NO_MPI_TIME[${d}]}" experiments/mpi/job.sh \
                  make -k -C experiments/mpi/aec run "${VARIABLES[@]}" D="${d}" "${FILES[@]}" \
                  THREADS="${AEC_THREADS}" TRIALS="${TRIALS}")
        echo "${DATASET}, ell = ${ELL}, d = ${d}, ALENEX27_AEC, ${AEC_THREADS} threads: job ${aec}"
        id=$(submit --dependency="afterany:${collect}:${aec}" -t 00:10:00 experiments/mpi/job.sh \
                 make -C experiments/mpi/aec merge "${VARIABLES[@]}" D="${d}")
        echo "${DATASET}, ell = ${ELL}, d = ${d}, merged with ALENEX27_AEC: job ${id}"
    fi
done
