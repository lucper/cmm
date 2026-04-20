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

    edges = [(u, v) for u in G for v in G[u] if u < v]

    for u, v in edges:
        patterns_u = set()
        patterns_v = set()

        for H in Hs:
            for i in range(len(u) - ell + 1):
                pattern = "".join(['*' if j in H else u[i+j] for j in range(ell)])
                patterns_u.add(pattern)
            for i in range(len(v) - ell + 1):
                pattern = "".join(['*' if j in H else v[i+j] for j in range(ell)])
                patterns_v.add(pattern)

        seen_in_edge = set()
        for x in patterns_u:
            for y in patterns_v:
                # Use a sorted tuple to treat {x, y} and {y, x} as the same key
                pair = tuple(sorted((x, y)))
                if pair not in seen_in_edge:
                    C[pair] += 1
                    seen_in_edge.add(pair)

    if not C:
        return {}

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
