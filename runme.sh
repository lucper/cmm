#!/usr/bin/bash

set -euo pipefail

WDIR="experiments"
OUTDIR="output"

## 1) Make every script in the tree executable (a fresh unzip may drop +x).
echo "=== Making scripts executable ===" >&2
find "${WDIR}" -type f \( -name '*.sh' -o -name '*.py' \) -exec chmod +x {} +

## 2) Compile the tools (cmm_test, cmm_eval, cmm_dedup) into experiments/bin.
echo "=== Building binaries ===" >&2
make -C "${WDIR}/src"

## 3) Prepare data: download + clean real data, sample artificial instances,
##    and unpack the bundled SLIDER datasets (shipped in data/).
echo "=== Preparing datasets ===" >&2
(
    cd "${WDIR}/data"
    ./run.sh -O "9606:Homo_sapiens 3702:Arabidopsis_thaliana 4932:Saccharomyces_cerevisiae 511145:Escherichia_coli_K12_MG1655 7227:Drosophila_melanogaster" \
             -s 700 -S 42 -x 9606 \
             -e "05 10 15 20 25 30" -v "100 200 400 800 1600 3200"
    tar zxvf slider_data.tar.gz
)

## 4) Unpack SLIDER (the heuristic jar) for the accuracy experiments.
echo "=== Unpacking SLIDER ===" >&2
(
    cd "${WDIR}/rq1-2_accuracy"
    if [ ! -f 000-slider/SliderLight/dist/SliderLight.jar ]; then
        unzip -o 000-slider/SliderLight.zip -d 000-slider
    fi
)

## ---------------------------------------------------------------------------
## RQ1 and RQ2 -- accuracy
## ---------------------------------------------------------------------------
echo "=== Running RQ1-RQ2 (accuracy) ===" >&2

P="64"          # threads for the exact miner
Ds="5"          # wildcard counts to sweep (single value for the reviewer portion)
STRING_DATA="511145 4932"
SLIDER_DATA="yeast_ii"

(
    cd "${WDIR}/rq1-2_accuracy"
    for d in ${Ds}
    do
        for dataset in ${STRING_DATA}
        do
            ./run.sh -D "${dataset}" -r 10 -t 10 -l 8 -d "${d}" -k 10000 -K 10000 -H 8 -c 70 -p "${P}" \
                     -f "../data/real_instances/${dataset}.sequences.s700.clean.fa" \
                     -i "../data/real_instances/${dataset}.links.s700.clean.int"
        done

        for dataset in ${SLIDER_DATA}
        do
            ./run.sh -D "${dataset}" -r 10 -t 10 -l 8 -d "${d}" -k 10000 -K 10000 -H 8 -c 70 -p "${P}" \
                     -f "../data/slider_data/${dataset}.clean.fa" \
                     -i "../data/slider_data/${dataset}.clean.int"
        done
    done
)

## ---------------------------------------------------------------------------
## RQ3, RQ4, RQ5 -- efficiency
## ---------------------------------------------------------------------------
echo "=== Running RQ3-RQ5 (efficiency) ===" >&2
(
    cd "${WDIR}/rq3-5_efficiency"
    ./run.sh -e "05 10 15 20 25 30" -v "100 200 400 800 1600 3200" -x 9606 \
             -p 128 -r 1 -k 10000 -l "5 8" -E 05 -L 8 -d 2
)

## ---------------------------------------------------------------------------
## Aggregate results into a single output/ tree.
## ---------------------------------------------------------------------------
echo "=== Aggregating results into ${OUTDIR}/ ===" >&2

ACC_OUT="${OUTDIR}/rq1-2_accuracy"
EFF_OUT="${OUTDIR}/rq3-5_efficiency"
mkdir -p "${ACC_OUT}" "${EFF_OUT}"

## Datasets that were run above (order defines the subfigure grid).
ALL_DATASETS="${STRING_DATA} ${SLIDER_DATA}"

## Auxiliary function for RQ1-2: assemble one multi-panel PDF.
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

## Auxiliary function for RQ3-5: Assemble the efficiency plots into one document.
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

        # rq3 and rq4: one figure page each, the two ell groups stacked.
        for rq in rq3 rq4; do
            for g in ell5_d0_1 ell5_d2_4 ell8_d0_3 ell8_d4_7; do
                _eff_block "${rq}_${g}"
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

## Make figures for RQ1-2.
build_figure_pdf "${WDIR}/rq1-2_accuracy/035-eval-density" \
                 "density.@.l8d5.h8.dedup" "density_legend" "${ACC_OUT}/density_all.pdf"
build_figure_pdf "${WDIR}/rq1-2_accuracy/030-score-curve" \
                 "curve.@.l8d5.dedup" "curve_legend" "${ACC_OUT}/curve_all.pdf"

cov_summary="${ACC_OUT}/coverage_summary.tsv"
printf 'dataset\tmethod\treduced_solution\tfull_solution\tcoverage\n' > "${cov_summary}"
cov_dir="${WDIR}/rq1-2_accuracy/035-eval-coverage"
if [ -d "${cov_dir}" ]; then
    find "${cov_dir}" -name '*.dat' | sort | while read -r dat; do
        fname=$(basename "${dat}")
        dataset=${fname%%.*}
        method=$(echo "${fname}" | cut -d. -f2)
        awk -v ds="${dataset}" -v m="${method}" \
            '{print ds "\t" m "\t" $1 "\t" $2 "\t" $3}' "${dat}" >> "${cov_summary}"
    done
fi

perf_summary="${ACC_OUT}/performance_summary.tsv"
printf 'dataset\ttime\tpeak_ram_mb_per_thread\ttotal_peak_ram_mb\n' > "${perf_summary}"
perf_dir="${WDIR}/rq1-2_accuracy/010-soln-cmm"
if [ -d "${perf_dir}" ]; then
    find "${perf_dir}" -mindepth 2 -type f -name '*.perf' | sort | while read -r perf_file; do
        # Extract the directory name directly containing the .perf file
        dataset=$(basename "$(dirname "${perf_file}")")

        awk -v ds="${dataset}" '
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

                printf "%s\t%02d:%02d:%02d\t%.2f\t%.2f\n", ds, hh, mm, ss, ram_mb_per_thread, ram_mb
            }
        ' "${perf_file}" >> "${perf_summary}"
    done
fi

## Make figures for RQ3-5.
if [ -d "${WDIR}/rq3-5_efficiency/plots" ]; then
    build_rq3-5_pdf "${WDIR}/rq3-5_efficiency/plots" "${EFF_OUT}/efficiency_all.pdf"
fi

echo "=== Done ===" >&2
echo "Aggregated results:" >&2
echo "  Accuracy   : ${ACC_OUT}/  (density_all.pdf, curve_all.pdf, coverage_summary.tsv)" >&2
echo "  Efficiency : ${EFF_OUT}/  (efficiency_all.pdf)" >&2
