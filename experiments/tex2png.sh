#!/usr/bin/bash
#
# Renders figures written by the gnuplot scripts of the experiments (cairolatex: a .tex with the
# text, and a .pdf with the graphics) to <figure>.png next to each .tex, for previewing. The
# figures themselves are not changed.
#
# Needs pdflatex and pdftoppm (poppler-utils).

set -euo pipefail

export LC_ALL=C

usage() {
    cat >&2 <<EOF
Usage: $0 [-r <dpi>] <figure>.tex ...

Options:
  -r <int>   resolution of the PNG, in dots per inch (default: 150)
  -h         show this help and exit
EOF
}

DPI=150
while getopts "r:h" opt
do
    case ${opt} in
        r) DPI=${OPTARG} ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done
shift $((OPTIND - 1))
if [ $# -eq 0 ]
then
    usage
    exit 1
fi

TMP=$(mktemp -d)
trap 'rm -rf "${TMP}"' EXIT

status=0
for tex in "$@"
do
    dir=$(cd "$(dirname "${tex}")" && pwd)
    name=$(basename "${tex}" .tex)

    # One page the size of the figure: the figure is boxed, and the page set to the box.
    printf '%s\n' \
        '\documentclass{article}' \
        '\usepackage{graphicx}' \
        '\usepackage{color}' \
        '\newsavebox\panel' \
        '\begin{document}' \
        "\\sbox\\panel{\\input{${name}}}" \
        '\pdfpagewidth=\wd\panel' \
        '\pdfpageheight=\dimexpr\ht\panel+\dp\panel\relax' \
        '\hoffset=-1in \voffset=-1in' \
        '\shipout\box\panel' \
        '\end{document}' > "${TMP}/page.tex"

    if ( cd "${TMP}" && TEXINPUTS="${dir}:" pdflatex -interaction=nonstopmode page.tex > page.out 2>&1 )
    then
        pdftoppm -png -r "${DPI}" -singlefile "${TMP}/page.pdf" "${dir}/${name}"
        echo "${dir}/${name}.png"
    else
        echo "Error: could not render ${tex}:" >&2
        grep -E '^!|not found' "${TMP}/page.log" | head -5 >&2 || true
        status=1
    fi
done
exit ${status}
