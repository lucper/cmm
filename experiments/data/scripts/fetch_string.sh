#!/usr/bin/env bash
#
# For each organism given on the command line (NCBI taxon ID + name), this
# script:
#   1. downloads the protein sequences (.fa.gz) and links (.links.txt.gz),
#   2. decompresses both,
#   3. for the given score cutoff, writes:
#        a) an unweighted undirected edge list (.int):
#        b) a matching FASTA (.fa) containing only proteins that still have
#           at least one edge at this cutoff (isolated nodes removed).
#
# Everything is written into <outdir>/<taxid>_<name>/ alongside the originals.
#
# STRING combined scores are scaled 0-1000:
#   150 = low, 400 = medium, 700 = high, 900 = highest confidence.
#   See https://string-db.org/help/scores/ and https://string-db.org/help/faq/
#
# Usage:
#   ./fetch_string.sh <threshold> <taxid:name> [<taxid:name> ...]
#
# Examples:
#   ./fetch_string.sh 700 9606:Homo_sapiens
#   ./fetch_string.sh 700 9606:Homo_sapiens 4932:Saccharomyces_cerevisiae
#   OUTDIR=string_data ./fetch_string.sh 400 9606:Homo_sapiens

set -euo pipefail

VERSION="v12.0"
BASE="https://stringdb-static.org/download"
OUTDIR="${OUTDIR:-string_data}"

usage() {
    echo "Usage: $0 <threshold> <taxid:name> [<taxid:name> ...]" >&2
    echo "  threshold  : STRING combined_score cutoff (0-1000), e.g. 150, 400, 700, 900" >&2
    echo "  taxid:name : NCBI taxon ID and organism name, e.g. 9606:Homo_sapiens" >&2
    echo "" >&2
    echo "Example: $0 700 9606:Homo_sapiens 4932:Saccharomyces_cerevisiae" >&2
    exit 1
}

if [ "$#" -lt 2 ]; then
    usage
fi

THRESHOLD="$1"
# drop first arg so that, from here, $@ contains the organisms only
shift

# Validate THRESHOLD
case "${THRESHOLD}" in
    ''|*[!0-9]*) echo "Error: threshold must be an integer in [0,1000]." >&2; exit 1 ;;
esac
[ "${THRESHOLD}" -le 1000 ] || { echo "Error: threshold must be <= 1000." >&2; exit 1; }

# Validate taxid:name pairs
ORGANISMS=("$@")
for entry in "${ORGANISMS[@]}"; do
    case "${entry}" in
        *:*) ;;
        *) echo "Error: organism must be given as 'taxid:name', got: ${entry}" >&2; exit 1 ;;
    esac
    taxid="${entry%%:*}"
    name="${entry#*:}"
    case "${taxid}" in
        ''|*[!0-9]*) echo "Error: taxid must be numeric, got: ${taxid}" >&2; exit 1 ;;
    esac
    [ -n "${name}" ] || { echo "Error: empty organism name in: ${entry}" >&2; exit 1; }
done

# Downloader: prefer wget, fall back to curl.
if command -v wget >/dev/null 2>&1; then
    fetch() { wget -q --show-progress -c -O "$1" "$2"; }
elif command -v curl >/dev/null 2>&1; then
    fetch() { curl -fL -C - -o "$1" "$2"; }
else
    echo "Error: need either wget or curl installed." >&2
    exit 1
fi

mkdir -p "${OUTDIR}"

for entry in "${ORGANISMS[@]}"; do
    taxid="${entry%%:*}"
    name="${entry#*:}"

    echo "=== ${name} (taxid ${taxid}) ===" >&2

    org_dir="${OUTDIR}/${taxid}_${name}"
    mkdir -p "${org_dir}"

    seq_gz="${org_dir}/${taxid}.protein.sequences.${VERSION}.fa.gz"
    links_gz="${org_dir}/${taxid}.protein.links.${VERSION}.txt.gz"
    seq_fa="${org_dir}/${taxid}.protein.sequences.${VERSION}.fa"
    links_txt="${org_dir}/${taxid}.protein.links.${VERSION}.txt"

    seq_url="${BASE}/protein.sequences.${VERSION}/${taxid}.protein.sequences.${VERSION}.fa.gz"
    links_url="${BASE}/protein.links.${VERSION}/${taxid}.protein.links.${VERSION}.txt.gz"

    # 1. Download (resume / skip if already present).
    [ -s "${seq_gz}" ]   || { echo "  downloading sequences ..." >&2; fetch "${seq_gz}"   "${seq_url}"; }
    [ -s "${links_gz}" ] || { echo "  downloading links ..."     >&2; fetch "${links_gz}" "${links_url}"; }

    # 2. Decompress (keep the .gz; -k, skip if already extracted).
    [ -s "${seq_fa}" ]    || { echo "  extracting sequences ..." >&2; gunzip -kf "${seq_gz}"; }
    [ -s "${links_txt}" ] || { echo "  extracting links ..."     >&2; gunzip -kf "${links_gz}"; }

    # 3. Build the unweighted undirected edge list AND a matching FASTA
    # containing only the proteins that still have at least one edge.
    int_out="${org_dir}/${taxid}.links.s${THRESHOLD}.int"
    fa_out="${org_dir}/${taxid}.sequences.s${THRESHOLD}.fa"

    if [ -s "${int_out}" ] && [ -s "${fa_out}" ]; then
        echo "  cutoff ${THRESHOLD}: already done, skipping" >&2
        continue
    fi

    echo "  cutoff ${THRESHOLD}: building edge list ..." >&2
    awk -v min="${THRESHOLD}" '
        NR == 1  { next }                       # header
        NF < 3   { next }                       # malformed
        $3 < min { next }                       # below cutoff
        {
            u = $1; v = $2;
            if (u > v) { t = u; u = v; v = t; } # canonical orientation
            key = u SUBSEP v;
            if (key in seen) next;              # collapse A B / B A
            seen[key] = 1;
            print u " " v;
        }
    ' "${links_txt}" > "${int_out}"

    echo "  cutoff ${THRESHOLD}: filtering FASTA to non-isolated nodes ..." >&2
    # First pass: collect the set of node IDs that appear in the edge list.
    # Second pass: stream the FASTA, keeping a record only if its header ID
    # is in that set.
    awk '
        FNR == NR { keep[$1] = 1; keep[$2] = 1; next }   # edge-list pass
        /^>/ {
            id = substr($1, 2);                          # strip ">"
            emit = (id in keep);
            if (emit) print;
            next;
        }
        { if (emit) print; }                             # sequence lines
    ' "${int_out}" "${seq_fa}" > "${fa_out}"

    n_edges=$(wc -l < "${int_out}")
    n_nodes=$(grep -c '^>' "${fa_out}" || true)
    echo "    cutoff ${THRESHOLD}: ${n_nodes} nodes, ${n_edges} edges" >&2
done

echo "" >&2
echo "Done. Each organism has its own directory under ./${OUTDIR}/ :" >&2
echo "  ${OUTDIR}/<taxid>_<name>/" >&2
echo "    {taxid}.sequences.s${THRESHOLD}.fa   (nodes, isolated proteins removed)" >&2
echo "    {taxid}.links.s${THRESHOLD}.int      (unweighted undirected edges)" >&2
echo "Pair the .fa and .int sharing the same sNNN suffix." >&2
