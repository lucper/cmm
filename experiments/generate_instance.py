#!/usr/bin/python3

from collections import defaultdict
import random
import sys
import string

def generate_graph(n, m, min_label_len, max_label_len, alphabet):
    ## Checks
    if min_label_len > max_label_len:
        raise ValueError(f'min_label_len ({min_label_len}) must be <= max_label_len ({max_label_len})')

    max_m = (n * (n - 1)) // 2
    if m > max_m:
        raise ValueError(f'Too many edges: max is {max_m} for {n} nodes.')

    min_m = (n + 1) // 2
    if m < min_m:
        raise ValueError(f"Too few edges: need at least {min_m} to avoid isolated vertices")

    max_possible_labels = sum(len(alphabet)**l for l in range(min_label_len, max_label_len + 1))
    if n > max_possible_labels:
        raise ValueError(f"Too many nodes: max is {max_possible_labels} for given alphabet and label lengths")
    ##

    ## Construct
    labels = set()
    while len(labels) < n:
        label_len = random.randint(min_label_len, max_label_len)
        labels.add("".join(random.choices(alphabet, k=label_len)))
    labels = sorted(list(labels))

    # Claude, Sonnet 4.6
    # Generate a random perfect (or near-perfect) matching to ensure no isolated vertices
    shuffled = labels[:]
    random.shuffle(shuffled)
    edges = set()
    for k in range(0, n - 1, 2):
        u, v = sorted([shuffled[k], shuffled[k + 1]])
        edges.add((u, v))
    # If num_nodes is odd, the last node is still isolated — connect it to a random neighbor
    if n % 2 == 1:
        u, v = sorted([shuffled[-1], random.choice(shuffled[:-1])])
        edges.add((u, v))
    # Fill up to m edges
    while len(edges) < m:
        u, v = sorted(random.sample(labels, 2))
        edges.add((u, v))

    G = {label: [] for label in labels}
    for u, v in edges:
        G[u].append(v)
        G[v].append(u)
    ##

    return G

if __name__ == '__main__':
    if len(sys.argv) != 9:
        raise SystemExit(f'Usage: {sys.argv[0]} <number of nodes> <number of edges> <minimum label length> <maximum label length> <alphabet size> <seed> <output_nodes_file> <output_edges_file>')
    n, m, min_label, max_label, alphabet_size, seed, nodes_output, edges_output = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), sys.argv[7], sys.argv[8]

    ascii_chars = string.ascii_lowercase + string.ascii_uppercase + string.digits
    if alphabet_size > len(ascii_chars):
        raise SystemExit(f'alphabet size ({alphabet_size}) must be at most {len(ascii_chars)} (ASCII lowercase, uppercase, and digits)')
    alphabet = ascii_chars[:alphabet_size]

    random.seed(seed)

    G = generate_graph(n, m, min_label, max_label, alphabet)
    V_ids = {v: i+1 for i, v in enumerate(G.keys())}

    with open(nodes_output, 'w') as fout:
        for i, node in enumerate(G.keys()):
            print(f'>{V_ids[node]}', file=fout)
            print(node, file=fout)

    with open(edges_output, 'w') as fout:
        for v in G:
            for u in G[v]:
                if v < u:
                    print(f'{V_ids[v]} {V_ids[u]}', file=fout)
