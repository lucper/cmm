#!/usr/bin/bash
#
# Runs the main experiment of the MPI version end to end. On a login node:
#   experiments/mpi/runme.sh                    # choose and edit what to run, then submit
#   experiments/mpi/runme.sh --yes              # submit the preset main as it is, without asking
#   experiments/mpi/runme.sh --dry-run          # print the sbatch commands, and submit nothing
# It submits the build, and one job per run, so that as many run at once as Slurm has nodes for: one per number
# of wildcards, of nodes, layout and trial, which wait for the build; for each number of wildcards one job that
# collects the results, which waits for those; and the code of the tag ALENEX27_AEC, merged with them: see README.md.
# A run that has finished gets no job, so submitting again runs only the others, and adding a layout runs only that
# layout. Turn the build off while jobs of this clone are running, as a build replaces their programs.
#
# The presets, of which the screen shows the parameters:
# - main: the layout 1x128 on 1, 2, 4 and 8 nodes, with the code of the tag ALENEX27_AEC (025-aec/).

set -euo pipefail

# The parameters, which a preset sets, and the screen edits.
declare -A JOBS AEC_TIME
PRESETS=(main)

# The layout 1x128 on 1, 2, 4 and 8 nodes, with ALENEX27_AEC.
preset_main() {
    PRESET=main
    DESCRIPTION="1x128 on 1, 2, 4 and 8 nodes, with ALENEX27_AEC"
    BUILD=1
    DATASET=human_ii
    FASTA=experiments/data/slider_data/human_ii.clean.fa    # from the root of the repository
    INTS=experiments/data/slider_data/human_ii.clean.int
    ELL=8
    SUPPORT=x2
    K=10000
    LAYOUTS="1x128"              # processes per node x threads per process; others: 8x16 2x64
    TRIALS=1
    AEC=1                        # also run the code of the tag ALENEX27_AEC, and merge it with the results (025-aec/)
    AEC_THREADS=128
    WILDCARDS="5"                # the numbers of wildcards to run, of those of JOBS
    # For each number of wildcards, the numbers of nodes, with the time limit of each run there: about twice what
    # the slowest layout, 1x128, needs, or more. At d = 5, it took 19, 11.5, 7.5, and 5.7 minutes on 1, 2, 4, and
    # 8 nodes. At d = 3, it took 11.3 hours on 1 node, 35 times as long; if the others scale as at d = 5, they need
    # about 7.4, 4.8, and 3.6 hours (a job of the 3 layouts on 8 nodes took 10% longer than such an estimate).
    JOBS=([5]="1:01:00:00 2:00:40:00 4:00:30:00 8:00:20:00"
          [3]="1:24:00:00 2:15:00:00 4:10:00:00 8:08:00:00")
    # For each of them, the time limit of a run of ALENEX27_AEC, on 1 node.
    AEC_TIME=([5]="01:00:00" [3]="24:00:00")
    EDITED=0
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

# The .perf file of a run: of a layout of d on nodes nodes, or of ALENEX27_AEC, in a trial.
mpi_perf() {
    echo "experiments/mpi/010-measure/${DATASET}/${DATASET}.layouts.l${ELL}d$1.${SUPPORT}.k${K}.n$2.$3.trial$4.perf"
}
aec_perf() {
    echo "experiments/mpi/025-aec/${DATASET}/${DATASET}.layouts.l${ELL}d$1.${SUPPORT}.k${K}.n1.aec${AEC_THREADS}.trial$2.perf"
}

plural() { [ "$1" = 1 ] && echo "$1 $2" || echo "$1 $2s"; }

# What a line of the plan says of its runs, of which todo of all have no results yet: the jobs it submits, and
# how many runs are done, which get none.
jobs_note() {
    local todo=$1 all=$2
    if [ "${todo}" = 0 ]; then echo "${GREEN}all $(plural "${all}" run) done${RESET}"
    elif [ "${todo}" = "${all}" ]; then plural "${todo}" job
    else echo "$(plural "${todo}" job), $((all - todo)) of $(plural "${all}" run) done"
    fi
}

# Prints the jobs that would be submitted, a line per group of runs, and sets NUM_JOBS and NUM_DONE.
plan() {
    NUM_JOBS=0 NUM_DONE=0
    local d job nodes time layout trial todo all
    if [ "${BUILD}" = 1 ]; then
        printf '  %-40s %-12s %s\n' "build" "" "build_snellius.sh"
        NUM_JOBS=$((NUM_JOBS + 1))
    fi
    local aec_todo=0
    if [ "${AEC}" = 1 ]; then
        for d in ${WILDCARDS}; do
            for trial in $(seq 1 "${TRIALS}"); do [ -e "$(aec_perf "${d}" "${trial}")" ] || aec_todo=$((aec_todo + 1)); done
        done
        if [ "${aec_todo}" != 0 ]; then
            printf '  %-40s %-12s %s\n' "build of ALENEX27_AEC" "00:30:00" "1 job"
            NUM_JOBS=$((NUM_JOBS + 1))
        fi
    fi
    for d in ${WILDCARDS}; do
        for job in ${JOBS[${d}]}; do
            nodes=${job%%:*} time=${job#*:}
            todo=0 all=0
            for trial in $(seq 1 "${TRIALS}"); do
                for layout in ${LAYOUTS}; do
                    all=$((all + 1)); [ -e "$(mpi_perf "${d}" "${nodes}" "${layout}" "${trial}")" ] || todo=$((todo + 1))
                done
            done
            printf '  %-40s %-12s %s\n' "       $(plural "${nodes}" node), ${LAYOUTS// /, }" "${time}" "$(jobs_note "${todo}" "${all}")"
            NUM_JOBS=$((NUM_JOBS + todo)) NUM_DONE=$((NUM_DONE + all - todo))
        done
        printf '  %-40s %-12s %s\n' "       collect the results (020-concat)" "00:10:00" "1 job, after those"
        NUM_JOBS=$((NUM_JOBS + 1))
        if [ "${AEC}" = 1 ]; then
            todo=0 all=0
            for trial in $(seq 1 "${TRIALS}"); do
                all=$((all + 1)); [ -e "$(aec_perf "${d}" "${trial}")" ] || todo=$((todo + 1))
            done
            printf '  %-40s %-12s %s\n' "       ALENEX27_AEC, ${AEC_THREADS} threads" "${AEC_TIME[${d}]}" "$(jobs_note "${todo}" "${all}")"
            printf '  %-40s %-12s %s\n' "       merge with ALENEX27_AEC (025-aec/)" "00:10:00" "1 job, after those"
            NUM_JOBS=$((NUM_JOBS + todo + 1)) NUM_DONE=$((NUM_DONE + all - todo))
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
    echo "${BOLD}The main experiment of the MPI version${RESET}  ${DIM}experiments/mpi/runme.sh${RESET}"
    echo
    field p "preset" "${BOLD}${PRESET}${RESET}${edited}  ${DIM}${DESCRIPTION}${RESET}"
    echo
    field 1 "build first" "$(yes_no "${BUILD}")"
    field 2 "name of the results" "${DATASET}"
    field 3 "sequences" "${FASTA}$(exists_note "${FASTA}")"
    field 4 "interactions" "${INTS}$(exists_note "${INTS}")"
    field 5 "motif length" "${ELL}"
    field 6 "wildcards" "${WILDCARDS}"
    field 7 "nodes:time per run" "${jobs%; }"
    field 8 "support function" "${SUPPORT}"
    field 9 "k" "${K}"
    field a "layouts" "${LAYOUTS}"
    field b "trials" "${TRIALS}"
    field c "ALENEX27_AEC" "$(yes_no "${AEC}")$([ "${AEC}" = 1 ] && echo ", ${AEC_THREADS} threads")"
    echo
    local lines
    lines=$(plan)
    plan > /dev/null
    echo "${BOLD}$(plural "${NUM_JOBS}" job)${RESET}, one per run, which Slurm runs as soon as it has the nodes:"
    echo "${lines}"
    if [ "${NUM_DONE}" != 0 ]; then
        echo "${DIM}A run that is done gets no job; to run it again, change the name of the results (2).${RESET}"
    fi
    local running
    running=$(running_jobs)
    if [ "${BUILD}" = 1 ] && [ "${running}" != 0 ]; then
        echo
        echo "${YELLOW}$(plural "${running}" job) of this clone queued or running: the build would replace their programs.${RESET}"
    fi
    [ -n "${MESSAGE:-}" ] && { echo; echo "${MESSAGE}"; MESSAGE=""; }
    echo
    echo "${DIM}p${RESET} preset  ${DIM}1-9 a-c${RESET} change  ${DIM}⏎${RESET} submit  ${DIM}n${RESET} dry run  ${DIM}q${RESET} quit"
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

# Reads the jobs and the time limit of ALENEX27_AEC of each number of wildcards.
edit_jobs() {
    local d value
    for d in ${WILDCARDS}; do
        value=${JOBS[${d}]:-1:01:00:00}
        read -r -e -i "${value}" -p "  d = ${d}, nodes:time limit of each run there: " value || return 0
        value=$(echo "${value}" | xargs)
        if [[ "${value}" =~ ^([0-9]+:[0-9:-]+)( [0-9]+:[0-9:-]+)*$ ]]; then
            [ "${value}" != "${JOBS[${d}]:-}" ] && EDITED=1
            JOBS[${d}]=${value}
        else
            MESSAGE="${RED}Not a list of nodes:time limit: '${value}'${RESET}"
            return 0
        fi
        value=${AEC_TIME[${d}]:-01:00:00}
        read -r -e -i "${value}" -p "  d = ${d}, time limit of ALENEX27_AEC: " value || return 0
        value=$(echo "${value}" | xargs)
        if [[ "${value}" =~ ^[0-9:-]+$ ]]; then
            [ "${value}" != "${AEC_TIME[${d}]:-}" ] && EDITED=1
            AEC_TIME[${d}]=${value}
        else
            MESSAGE="${RED}Not a time limit: '${value}'${RESET}"
        fi
    done
}

# Every number of wildcards needs its jobs: asks for those of the new ones.
check_wildcards() {
    local d
    for d in ${WILDCARDS}; do
        [ -n "${JOBS[${d}]:-}" ] && [ -n "${AEC_TIME[${d}]:-}" ] || { edit_jobs; return; }
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
            c) toggle AEC; [ "${AEC}" = 1 ] && edit AEC_THREADS "threads of ALENEX27_AEC" "${INT}" ;;
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
    [ -n "${JOBS[${d}]:-}" ] && [ -n "${AEC_TIME[${d}]:-}" ] || { echo "No jobs for d = ${d}" >&2; exit 2; }
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

# The build of ALENEX27_AEC, unless all its runs are done.
aec_todo=0
if [ "${AEC}" = 1 ]; then
    for d in ${WILDCARDS}; do
        for trial in $(seq 1 "${TRIALS}"); do [ -e "$(aec_perf "${d}" "${trial}")" ] || aec_todo=$((aec_todo + 1)); done
    done
fi
aec_build=""
if [ "${aec_todo}" != 0 ]
then
    aec_build=$(submit --ntasks=1 --cpus-per-task=16 -t 00:30:00 experiments/mpi/job.sh \
                    make -C experiments/mpi/025-aec build)
    echo "build of ALENEX27_AEC: job ${aec_build}"
fi

for d in ${WILDCARDS}
do
    FILES=(FASTA="${PWD}/${FASTA}" INTS="${PWD}/${INTS}")
    ids=""    # the runs of d, for which the collection waits
    skipped=0

    # One job per run, so that they all can run at once.
    for job in ${JOBS[${d}]}
    do
        nodes=${job%%:*} time=${job#*:}
        for trial in $(seq 1 "${TRIALS}")
        do
            for layout in ${LAYOUTS}
            do
                if [ -e "$(mpi_perf "${d}" "${nodes}" "${layout}" "${trial}")" ]; then skipped=$((skipped + 1)); continue; fi
                id=$(submit "${dependency[@]}" --nodes="${nodes}" --exclusive -t "${time}" experiments/mpi/job.sh \
                         make -k -C experiments/mpi/010-measure all "${VARIABLES[@]}" D="${d}" "${FILES[@]}" \
                         NODES="${nodes}" LAYOUTS="${layout}" TRIALS="${TRIALS}" TRIAL="${trial}")
                echo "${DATASET}, ell = ${ELL}, d = ${d}, ${nodes} node(s), ${layout}, trial ${trial}, ${time}: job ${id}"
                ids="${ids}:${id}"
            done
        done
    done
    [ "${skipped}" = 0 ] || echo "${DATASET}, ell = ${ELL}, d = ${d}: $(plural "${skipped}" run) done already, without a job"

    # afterany: also if a run failed or reached its limit, with the runs that finished.
    after=()
    [ -n "${ids}" ] && after=("--dependency=afterany${ids}")
    collect=$(submit "${after[@]}" -t 00:10:00 experiments/mpi/job.sh \
                  make -C experiments/mpi/020-concat all "${VARIABLES[@]}" D="${d}")
    echo "${DATASET}, ell = ${ELL}, d = ${d}, results: job ${collect}"

    # The code of the tag ALENEX27_AEC, apart from the above; the merge waits for both.
    if [ "${AEC}" = 1 ]
    then
        aec_ids=""
        for trial in $(seq 1 "${TRIALS}")
        do
            [ -e "$(aec_perf "${d}" "${trial}")" ] && continue
            aec=$(submit --dependency="afterok:${aec_build}" --nodes=1 --exclusive -t "${AEC_TIME[${d}]}" experiments/mpi/job.sh \
                      make -k -C experiments/mpi/025-aec run "${VARIABLES[@]}" D="${d}" "${FILES[@]}" \
                      THREADS="${AEC_THREADS}" TRIALS="${TRIALS}" TRIAL="${trial}")
            echo "${DATASET}, ell = ${ELL}, d = ${d}, ALENEX27_AEC, ${AEC_THREADS} threads, trial ${trial}: job ${aec}"
            aec_ids="${aec_ids}:${aec}"
        done
        id=$(submit --dependency="afterany:${collect}${aec_ids}" -t 00:10:00 experiments/mpi/job.sh \
                 make -C experiments/mpi/025-aec merge "${VARIABLES[@]}" D="${d}")
        echo "${DATASET}, ell = ${ELL}, d = ${d}, merged with ALENEX27_AEC: job ${id}"
    fi
done
