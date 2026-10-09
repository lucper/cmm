#!/usr/bin/bash
#
# Reproduces the experiments of the paper (RQ1-RQ5) and collects the results
# in output/full/ (or output/subset/ with -s). All parameters are set in
# sections 1-3 of this script:
#   1. Parameters shared by the subset and the full experiments
#   2. Subset of parameters used by the ALENEX 2027 AEC (-s)
#   3. Parameters of the full experiments of the paper (default)
#   4. Build
#   5. Data preparation
#   6. RQ1-RQ2 (accuracy)
#   7. RQ3-RQ5 (efficiency)
#   8. Aggregation of the results

set -euo pipefail

export LC_ALL=C

usage() {
    cat >&2 <<EOF
Usage: $0 [-s]

Options:
  -s  run the subset of the experiments evaluated by the ALENEX 2027 AEC
      instead of the full experiments of the paper
  -h  show this help and exit
EOF
}

SUBSET=0

while getopts "sh" opt
do
    case ${opt} in
        s) SUBSET=1 ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

## ---------------------------------------------------------------------------
## 1. Parameters shared by the subset and the full experiments
## ---------------------------------------------------------------------------
# Data preparation
ORGANISMS="9606:Homo_sapiens 3702:Arabidopsis_thaliana 4932:Saccharomyces_cerevisiae 511145:Escherichia_coli_K12_MG1655 7227:Drosophila_melanogaster"
THRESHOLD=700         # STRING combined score cutoff
SEED=42               # seed for sampling the artificial instances
INSTANCE_TAXID=9606   # organism the artificial instances are sampled from

# RQ1-RQ2
ELL=8                 # motif length
R=10                  # number of SLIDER runs
T=10                  # max convergence time of SLIDER (min)
K_SLIDER=10000        # SLIDER output size
K_TOP=10000           # kept top motif pairs after aggregation
H=8                   # similarity proximity
C=70                  # coverage cutoff

# RQ3-RQ5
K_CMM=10000           # number of top motif pairs kept by cmm_test
TRIALS=1              # trials per configuration
FIXED_ED=05           # edge density for RQ3 and RQ5
FIXED_V=200           # node count for RQ4 and RQ5
RQ5_D=2               # number of wildcards for RQ5

if [ ${SUBSET} -eq 1 ]
then
## ---------------------------------------------------------------------------
## 2. Subset of parameters used by the ALENEX 2027 AEC
## ---------------------------------------------------------------------------
    THREADS=64                           # threads for the exact algorithm
    DS="5"                               # RQ1-RQ2: numbers of wildcards
    STRING_DATA="511145 4932"            # RQ1-RQ2: datasets from STRING
    SLIDER_DATA="yeast_ii"               # RQ1-RQ2: datasets from the SLIDER evaluation
    DENSITIES="05 10 15 20 25"           # RQ3-RQ5: edge densities (%)
    NODE_COUNTS="100 200 400 800 1600"   # RQ3-RQ5: node counts
    ELLS="5"                             # RQ3-RQ4: motif lengths
    RQ5_ELL=5                            # RQ5: motif length
    OUTDIR="output/subset"               # where the results are collected
else
## ---------------------------------------------------------------------------
## 3. Parameters of the full experiments of the paper
## ---------------------------------------------------------------------------
    THREADS=128
    DS="3 5"
    STRING_DATA="9606 3702 4932 511145 7227"
    SLIDER_DATA="yeast_ii human_ii"
    DENSITIES="05 10 15 20 25 30"
    NODE_COUNTS="100 200 400 800 1600 3200"
    ELLS="5 8"
    RQ5_ELL=8
    OUTDIR="output/full"
fi

## ---------------------------------------------------------------------------
## 4. Build
## ---------------------------------------------------------------------------
WDIR="experiments"

# Make every script in the tree executable (a fresh unzip may drop +x).
echo "=== Making scripts executable ===" >&2
find "${WDIR}" -type f \( -name '*.sh' -o -name '*.py' \) -exec chmod +x {} +

# Compile the tools (cmm_test, cmm_eval, cmm_dedup) into experiments/bin.
echo "=== Building binaries ===" >&2
make -C "${WDIR}/src"

## ---------------------------------------------------------------------------
## 5. Data preparation
## ---------------------------------------------------------------------------
# Download and clean the real datasets, sample the artificial instances, and
# unpack the bundled datasets from the SLIDER evaluation.
echo "=== Preparing datasets ===" >&2
(
    cd "${WDIR}/data"
    ./run.sh -O "${ORGANISMS}" -s ${THRESHOLD} -S ${SEED} -x ${INSTANCE_TAXID} \
             -e "${DENSITIES}" -v "${NODE_COUNTS}"
    tar zxvf slider_data.tar.gz
)

# Unpack SLIDER (the heuristic jar) for the accuracy experiments.
echo "=== Unpacking SLIDER ===" >&2
(
    cd "${WDIR}/rq1-2_accuracy"
    if [ ! -f 000-slider/SliderLight/dist/SliderLight.jar ]; then
        unzip -o 000-slider/SliderLight.zip -d 000-slider
    fi
)

## ---------------------------------------------------------------------------
## 6. RQ1-RQ2 (accuracy)
## ---------------------------------------------------------------------------
echo "=== Running RQ1-RQ2 (accuracy) ===" >&2
(
    cd "${WDIR}/rq1-2_accuracy"
    for dataset in ${STRING_DATA}
    do
        ./run.sh -D "${dataset}" -r ${R} -t ${T} -l ${ELL} -d "${DS}" -k ${K_SLIDER} -K ${K_TOP} \
                 -H ${H} -c ${C} -p ${THREADS} \
                 -f "../data/real_instances/${dataset}.sequences.s${THRESHOLD}.clean.fa" \
                 -i "../data/real_instances/${dataset}.links.s${THRESHOLD}.clean.int"
    done

    for dataset in ${SLIDER_DATA}
    do
        ./run.sh -D "${dataset}" -r ${R} -t ${T} -l ${ELL} -d "${DS}" -k ${K_SLIDER} -K ${K_TOP} \
                 -H ${H} -c ${C} -p ${THREADS} \
                 -f "../data/slider_data/${dataset}.clean.fa" \
                 -i "../data/slider_data/${dataset}.clean.int"
    done
)

## ---------------------------------------------------------------------------
## 7. RQ3-RQ5 (efficiency)
## ---------------------------------------------------------------------------
echo "=== Running RQ3-RQ5 (efficiency) ===" >&2
(
    cd "${WDIR}/rq3-5_efficiency"
    ./run.sh -e "${DENSITIES}" -v "${NODE_COUNTS}" -x ${INSTANCE_TAXID} \
             -p ${THREADS} -r ${TRIALS} -k ${K_CMM} -l "${ELLS}" -E ${FIXED_ED} -V ${FIXED_V} -L ${RQ5_ELL} -d ${RQ5_D}
)

## ---------------------------------------------------------------------------
## 8. Aggregation of the results
## ---------------------------------------------------------------------------
echo "=== Aggregating results into ${OUTDIR}/ ===" >&2

mkdir -p "${OUTDIR}"

# Datasets that were run above (order defines the subfigure grid).
ALL_DATASETS="${STRING_DATA} ${SLIDER_DATA}"

# RQ1-RQ2: assembles one multi-panel PDF.
build_figure_pdf() {
    local sdir="$1" tmpl="$2" legend="$3" out="$4"
    local wrap; wrap=$(mktemp -d)

    # Bring the fragments + their graphics PDFs into the build dir.
    cp -f "${sdir}/"*.tex "${wrap}/" 2>/dev/null || true
    cp -f "${sdir}/"*.pdf "${wrap}/" 2>/dev/null || true

    {
        cat <<'HEAD'
\documentclass{article}
\usepackage[paperwidth=18cm,paperheight=26cm,margin=1cm]{geometry}
\pagestyle{empty}
\usepackage{graphicx}
\usepackage{color}
\usepackage{caption}
% Placeholder for a figure that was not generated. Uses only \framebox/\parbox
% (base LaTeX) so it needs no extra packages.
\newcommand{\missingfig}[1]{%
  \framebox[3.2cm]{\parbox[c][2.4cm][c]{3.0cm}{\centering\footnotesize\color{gray}#1\\[2pt](not generated)}}%
}
\newcommand{\subfig}[2]{%
  \begin{minipage}[t]{0.48\linewidth}
    \centering
    \IfFileExists{#1.tex}{\input{#1.tex}}{\missingfig{\texttt{\detokenize{#1}}}}
    \par\vspace{0.3em}
    \captionof{figure}{#2}
  \end{minipage}%
}
\begin{document}
\begin{minipage}{15cm}
  \centering
HEAD

        if [ -f "${wrap}/${legend}.tex" ]; then
            echo "  \\input{${legend}.tex}\\par\\vspace{0.4em}"
        fi

        local i=0
        for ds in ${ALL_DATASETS}; do
            local stem="${tmpl/@/${ds}}"
            printf '  \\subfig{%s}{\\texttt{\\detokenize{%s}}}%%\n' "${stem}" "${ds}"
            i=$((i + 1))
            if [ $((i % 2)) -eq 0 ]; then echo "  \\par\\vspace{0.6em}"; else echo "  \\hfill"; fi
        done

        cat <<'FOOT'
\end{minipage}
\end{document}
FOOT
    } > "${wrap}/figure.tex"

    ( cd "${wrap}" && pdflatex -interaction=nonstopmode figure.tex >figure.build.log 2>&1 ) || true

    if [ -f "${wrap}/figure.pdf" ]; then
        cp "${wrap}/figure.pdf" "${out}"
        rm -rf "${wrap}"
    else
        echo "  ERROR: could not build ${out}. LaTeX errors:" >&2
        grep -E '^!|Error|not found|Cannot|No file' "${wrap}/figure.build.log" | head -15 >&2 || true
        echo "  (full log and build dir kept at: ${wrap})" >&2
    fi
}

# RQ3-RQ5: assembles the efficiency plots into one PDF.
build_rq3-5_pdf() {
    local plots="$1" out="$2"
    local wrap; wrap=$(mktemp -d)
    cp -f "${plots}/"*.pdf "${wrap}/" 2>/dev/null || true

    # One legend+time+memory block for a given basename prefix.
    _eff_block() {
        local p="$1"
        cat <<BLOCK
\IfFileExists{${p}_legend.pdf}{\centerline{\includegraphics[width=0.6\linewidth]{${p}_legend.pdf}}\\\\[2pt]}{}%
\begin{minipage}[t]{0.49\linewidth}\centering
  \IfFileExists{${p}_time.pdf}{\includegraphics[width=\linewidth]{${p}_time.pdf}}{\effmissing{${p}_time}}%
\end{minipage}\hfill%
\begin{minipage}[t]{0.49\linewidth}\centering
  \IfFileExists{${p}_memory.pdf}{\includegraphics[width=\linewidth]{${p}_memory.pdf}}{\effmissing{${p}_memory}}%
\end{minipage}\par\vspace{0.2cm}
BLOCK
    }

    {
        cat <<'HEAD'
\documentclass{article}
\usepackage[paperwidth=18cm,paperheight=26cm,margin=1cm]{geometry}
\pagestyle{empty}
\usepackage{graphicx}
\usepackage{color}
\newcommand{\effmissing}[1]{\framebox[6cm]{\parbox[c][3.5cm][c]{5.5cm}{\centering\footnotesize\color{gray}\texttt{\detokenize{#1}}\\[2pt](not generated)}}}
\begin{document}\centering
HEAD

        # rq3 and rq4: one figure page each, with two d groups per ell (as in
        # d_groups of plot.py).
        for rq in rq3 rq4; do
            for ell in ${ELLS}; do
                local half=$((ell / 2))
                _eff_block "${rq}_ell${ell}_d0_$((half - 1))"
                _eff_block "${rq}_ell${ell}_d${half}_$((ell - 1))"
            done
            echo '\newpage'
        done
        # rq5: single thread-scaling block.
        _eff_block "rq5"

        echo '\end{document}'
    } > "${wrap}/efficiency.tex"

    ( cd "${wrap}" && pdflatex -interaction=nonstopmode efficiency.tex >efficiency.build.log 2>&1 ) || true

    if [ -f "${wrap}/efficiency.pdf" ]; then
        cp "${wrap}/efficiency.pdf" "${out}"
        rm -rf "${wrap}"
    else
        echo "  ERROR: could not build ${out}. LaTeX errors:" >&2
        grep -E '^!|Error|not found|Cannot|No file' "${wrap}/efficiency.build.log" | head -15 >&2 || true
        echo "  (full log and build dir kept at: ${wrap})" >&2
    fi
}

# RQ1-RQ2: figures, one PDF per number of wildcards.
for d in ${DS}; do
    build_figure_pdf "${WDIR}/rq1-2_accuracy/035-eval-density" \
                     "density.@.l${ELL}d${d}.h${H}.dedup" "density_legend" "${OUTDIR}/rq2_similarity_densities_l${ELL}d${d}.pdf"
    build_figure_pdf "${WDIR}/rq1-2_accuracy/030-score-curve" \
                     "curve.@.l${ELL}d${d}.dedup" "curve_legend" "${OUTDIR}/rq1_score_curves_l${ELL}d${d}.pdf"
done

# RQ1-RQ2: coverage table.
cov_summary="${OUTDIR}/rq2_coverage.tsv"
cov_dir="${WDIR}/rq1-2_accuracy/035-eval-coverage"
printf 'dataset\td\tmethod\treduced_solution\tfull_solution\tcoverage\n' > "${cov_summary}"
for ds in ${ALL_DATASETS}; do
    for d in ${DS}; do
        for method in seq_slider m_slider; do
            dat="${cov_dir}/${ds}/${ds}.${method}.l${ELL}d${d}.recall_c${C}.trials${R}.min${T}.k${K_TOP}.h${H}.dedup.dat"
            if [ -f "${dat}" ]; then
                awk -v ds="${ds}" -v d="${d}" -v m="${method}" \
                    '{print ds "\t" d "\t" m "\t" $1 "\t" $2 "\t" $3}' "${dat}" >> "${cov_summary}"
            else
                echo "  WARNING: missing ${dat}" >&2
            fi
        done
    done
done

# RQ1-RQ2: performance table of the exact algorithm.
perf_summary="${OUTDIR}/exact_algorithm_performance.tsv"
perf_dir="${WDIR}/rq1-2_accuracy/010-soln-cmm"
printf 'dataset\td\ttime\tpeak_ram_mb_per_thread\ttotal_peak_ram_mb\n' > "${perf_summary}"
for ds in ${ALL_DATASETS}; do
    for d in ${DS}; do
        perf_file="${perf_dir}/${ds}/${ds}.cmm.l${ELL}d${d}.k${K_TOP}.p${THREADS}.perf"
        if [ ! -f "${perf_file}" ]; then
            echo "  WARNING: missing ${perf_file}" >&2
            continue
        fi
        awk -v ds="${ds}" -v d="${d}" '
            NF >= 14 {
                # Convert milliseconds to seconds
                total_sec = int($13 / 1000)
                hh = int(total_sec / 3600)
                mm = int((total_sec % 3600) / 60)
                ss = total_sec % 60

                # Convert KB to MB
                ram_mb = $14 / 1024

                # Per thread
                ram_mb_per_thread = ram_mb / $12

                printf "%s\t%s\t%02d:%02d:%02d\t%.2f\t%.2f\n", ds, d, hh, mm, ss, ram_mb_per_thread, ram_mb
            }
        ' "${perf_file}" >> "${perf_summary}"
    done
done

# RQ3-RQ5: figures.
if [ -d "${WDIR}/rq3-5_efficiency/plots/p${THREADS}" ]; then
    build_rq3-5_pdf "${WDIR}/rq3-5_efficiency/plots/p${THREADS}" "${OUTDIR}/rq3-5_efficiency.pdf"
fi

echo "=== Done ===" >&2
echo "Aggregated results in ${OUTDIR}/ (see README.md)" >&2
