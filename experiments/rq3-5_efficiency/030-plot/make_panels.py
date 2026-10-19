#!/usr/bin/env python3
"""Split one TSV of 020-concat into the data files of its panels.

Usage: make_panels.py <rq> <tsv> <stem> [annotate]
  e.g. make_panels.py rq3 ../020-concat/rq3.9606.s700.S42.ed05.k10000.p128.tsv rq3.9606.s700.S42.ed05.k10000.p128 "8:4 8:7"

rq:        rq3, rq4, or rq5
annotate:  "ell:d ..." pairs whose time points are labelled with max_rank (RQ3, RQ4)
Input:     trial, then the 13 columns of cmm_test (see the README)
Output:    for RQ3 and RQ4, per motif length L, <stem>_ell<L>_d0_<L/2-1>.dat with the numbers of
           wildcards d < L/2 and <stem>_ell<L>_d<L/2>_<L-1>.dat with the others (L/2 rounded down), as
           in the plots of the paper; for RQ5, <stem>.dat. Their names are printed on stdout. Each
           holds one gnuplot data set per d (two blank lines apart), headed by its legend title, with
           one row per point:
           x <tab> N <tab> time_s <tab> mem_mb <tab> label
           x is V (RQ3), the edge density in percent (RQ4), or the number of threads (RQ5); time and
           memory are medians over the trials, memory per thread on RQ3 and RQ4; label is max_rank on
           annotated pairs, the speedup over the first point on RQ5, and "-" otherwise.
"""
import csv, sys
from statistics import median

if len(sys.argv) not in (4, 5):
    sys.exit("usage: make_panels.py <rq> <tsv> <stem> [annotate]")
rq, tsv, stem = sys.argv[1:4]
annotate = set(sys.argv[4].split()) if len(sys.argv) > 4 else set()

def x_of(row):
    return {"rq3": int(row[3]), "rq4": round(float(row[2]) * 100), "rq5": int(row[10])}[rq]

# (ell, d) -> x -> rows of its trials
points = {}
with open(tsv, newline="") as fh:
    for row in csv.reader(fh, delimiter="\t"):
        if row:
            points.setdefault((int(row[5]), int(row[6])), {}).setdefault(x_of(row), []).append(row)

# The lower and upper halves of the numbers of wildcards of ell, as in the plots of the paper.
def group(ell, d):
    if rq == "rq5":
        return stem
    half = ell // 2
    return f"{stem}_ell{ell}_d0_{half - 1}" if d < half else f"{stem}_ell{ell}_d{half}_{ell - 1}"

files = {}
for ell, d in sorted(points):
    name = group(ell, d)
    out = files.setdefault(name, [])
    if out:
        out.append("\n\n")
    # gnuplot reads backslash escapes in quoted strings of data files, hence the doubled backslashes.
    out.append(f'"$(\\\\ell={ell},\\\\ d={d})$"\n')
    base = None
    for x, rows in sorted(points[ell, d].items()):
        time_s = median(float(r[12]) for r in rows) / 1000
        mem_mb = median(float(r[13]) / (1 if rq == "rq5" else float(r[11])) for r in rows) / 1024
        base = base or time_s
        if rq == "rq5":
            label = f"{base / time_s:.1f}x"
        elif f"{ell}:{d}" in annotate:
            label = f"{int(rows[0][7]):,}"
        else:
            label = "-"
        out.append(f"{x}\t{rows[0][4]}\t{time_s:.15g}\t{mem_mb:.15g}\t{label}\n")

for name, out in files.items():
    with open(f"{name}.dat", "w") as fh:
        fh.writelines(out)
    print(f"{name}.dat")
