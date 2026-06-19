#!/usr/bin/env python3
"""Normalize one dataset's ranked score files into a tidy TSV on stdout.

Usage: transform.py <exact_token> <methods.csv> <file>...
  e.g. transform.py exact methods.csv score/ecoli/ecoli.*l8d3*dat > tidy/ecoli.tsv

Each input line:  COL1 COL2 score   (already sorted, best first)
methods.csv:      raw,"Proper"
Output (stdout):  rank <tab> proper_method <tab> normscore
"""
import csv, sys, os

exact_token = sys.argv[1]
methods_csv = sys.argv[2]
files       = sys.argv[3:]
if not files:
    sys.exit("usage: transform.py <exact_token> <methods.csv> <file>...")

mname = {}
with open(methods_csv, newline="") as fh:
    for row in csv.reader(fh):
        if row:
            mname[row[0]] = row[1]

# method -> list of scores (file is pre-sorted, keep order)
scores = {}
for f in files:
    method = os.path.basename(f).split(".")[1]
    scores[method] = [float(line.split()[2]) for line in open(f) if line.strip()]

if exact_token not in scores:
    sys.exit(f"no exact file ('{exact_token}') among inputs")

n = min(len(v) for v in scores.values())
scores = {m: v[:n] for m, v in scores.items()} # truncate to solution set with smallest items

opt = scores[exact_token][0] # top line of exact = optimum

out = sys.stdout
for method, vals in scores.items():
    proper = mname.get(method, method)
    for rank, s in enumerate(vals, 1):
        out.write(f"{rank}\t{proper}\t{s/opt:.5f}\n")
