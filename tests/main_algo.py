#!/usr/bin/python3

from itertools import product, combinations
from collections import defaultdict
import sys

def main_algo(G, ell, d):
    """Given a graph
        G = {'aba': ['abba', 'bba'],
             'bba': ['aba', 'abba'],
             'abba': ['bba', 'aba']}
    and a positive integer ell and non-negative integer d < ell,
    return the pairs of length-ell strings (x,y) with d wildcards 
    each such that the number of edges {u,v} with x ~< u and y ~< v
    is maximized. (Notation: w ~< z means that w is a substring of
    z with the wildcards positions matching any character.)
    """
    if d >= ell:
        raise ValueError(f'd ({d}) must be < ell ({ell})')
    Hs = list(combinations(range(ell), d))
    C = defaultdict(int)
    for H_x in Hs:
        for H_y in Hs:
            for v in G:
                V = {x for i in range(len(v)-ell+1)
                     for x in ["".join(['*' if j in H_x else v[i+j] for j in range(ell)])]}
                for u in G[v]:
                    if v < u:
                        U = {y for i in range(len(u)-ell+1)
                             for y in ["".join(['*' if j in H_y else u[i+j] for j in range(ell)])]}
                        for x, y in product(V, U):
                            C[x, y] += 1
    max_count = max(C.values())
    return {k: v for k, v in C.items() if v == max_count}

if __name__ == '__main__':
    if len(sys.argv) != 5:
        raise SystemExit(f'Usage: {sys.argv[0]} <ell> <d> <input_nodes_file> <input_edges_file>')
    ell, d, nodes_input, edges_input = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]

    V = defaultdict(str)
    with open(nodes_input, 'r') as fin:
        next(fin) # skip header
        for line in fin:
            i, label = line.split(';')
            V[int(i)] = label.strip()

    G = defaultdict(list)
    with open(edges_input, 'r') as fin:
        next(fin)
        for line in fin:
            u, v = line.split(';')
            G[V[int(u)]].append(V[int(v)])
            G[V[int(v)]].append(V[int(u)])

    counts = main_algo(G, ell, d)
    for motif, count in counts.items():
        x, y = motif
        print(x, y, count)
