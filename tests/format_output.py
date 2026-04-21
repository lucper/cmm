#!/usr/bin/python3

import fileinput
from pprint import pprint

output = []

for line in fileinput.input():
    X, Y, *counts = line.strip().split()
    if X <= Y:
        output.append((X,Y,"\t".join(counts)))
    else:
        output.append((Y,X,"\t".join(counts)))

output.sort()

for X, Y, counts in output:
    print(X, Y, counts, sep='\t')
