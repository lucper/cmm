#!/usr/bin/python3

import random
import sys
import string

AA_ALPHABET = list("ACDEFGHIKLMNPQRSTVWY")

def generate_graph(n, density, min_label_len, max_label_len, alphabet):
    if min_label_len > max_label_len:
        raise ValueError(f'min_label_len ({min_label_len}) must be <= max_label_len ({max_label_len})')

    if not (0.0 <= density <= 1.0):
        raise ValueError(f'density ({density}) must be in [0.0, 1.0]')

    max_possible_edges = (n * (n - 1)) // 2
    m = round(density * max_possible_edges)

    min_m = (n + 1) // 2
    if m < min_m:
        raise ValueError(
            f'Requested density {density} yields {m} edges, but at least {min_m} are needed '
            f'to avoid isolated vertices (try a higher density).'
        )

    max_possible_labels = sum(len(alphabet) ** l for l in range(min_label_len, max_label_len + 1))
    if n > max_possible_labels:
        raise ValueError(f'Too many nodes: max is {max_possible_labels} for given alphabet and label lengths')

    ## Construct labels
    labels = set()
    while len(labels) < n:
        label_len = random.randint(min_label_len, max_label_len)
        labels.add("".join(random.choices(alphabet, k=label_len)))
    labels = sorted(list(labels))

    # Start with a random perfect (or near-perfect) matching to ensure no isolated vertices.
    shuffled = labels[:]
    random.shuffle(shuffled)
    edges = set()
    for i in range(0, n - 1, 2):
        u, v = sorted([shuffled[i], shuffled[i + 1]])
        edges.add((u, v))
    # If n is odd, the last node is still isolated — connect it to a random neighbour.
    if n % 2 == 1:
        u, v = sorted([shuffled[-1], random.choice(shuffled[:-1])])
        edges.add((u, v))

    # Fill up to m edges by sampling pairs uniformly at random.
    # Uniform random sampling avoids any community or bipartite structure.
    if m > max_possible_edges:
        m = max_possible_edges

    while len(edges) < m:
        u, v = sorted(random.sample(labels, 2))
        edges.add((u, v))

    G = {label: [] for label in labels}
    for u, v in edges:
        G[u].append(v)
        G[v].append(u)

    return G


if __name__ == '__main__':
    if len(sys.argv) != 9:
        raise SystemExit(
            f'Usage: {sys.argv[0]} <number of nodes> <edge density [0,1]> '
            f'<minimum label length> <maximum label length> <alphabet size> '
            f'<seed> <output_nodes_file> <output_edges_file>'
        )

    n             = int(sys.argv[1])
    density       = float(sys.argv[2])
    min_label     = int(sys.argv[3])
    max_label     = int(sys.argv[4])
    alphabet_size = int(sys.argv[5])
    seed          = int(sys.argv[6])
    nodes_output  = sys.argv[7]
    edges_output  = sys.argv[8]

    if alphabet_size > len(AA_ALPHABET):
        raise SystemExit(
            f'alphabet size ({alphabet_size}) must be at most {len(AA_ALPHABET)} '
            f'(amino acid alphabet)'
        )
    alphabet = AA_ALPHABET[:alphabet_size]

    random.seed(seed)

    G = generate_graph(n, density, min_label, max_label, alphabet)
    V_ids = {v: i + 1 for i, v in enumerate(G.keys())}

    with open(nodes_output, 'w') as fout:
        for node in G.keys():
            print(f'>{V_ids[node]}', file=fout)
            print(node, file=fout)

    with open(edges_output, 'w') as fout:
        for v in G:
            for u in G[v]:
                if v < u:
                    print(f'{V_ids[v]} {V_ids[u]}', file=fout)
