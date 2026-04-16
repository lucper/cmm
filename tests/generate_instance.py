#!/usr/bin/python3

from collections import defaultdict
import random
import sys

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
    labels = list(labels)

    # Claude, Sonnet 4.6
    # Generate a random perfect (or near-perfect) matching to ensure no isolated vertices
    shuffled = labels[:]
    random.shuffle(shuffled)
    edges = set()
    for k in range(0, n - 1, 2):
        u, v = shuffled[k], shuffled[k + 1]
        if u > v:
            u, v = v, u
        edges.add((u, v))
    # If num_nodes is odd, the last node is still isolated — connect it to a random neighbor
    if n % 2 == 1:
        u = shuffled[-1]
        v = random.choice(shuffled[:-1])
        if u > v:
            u, v = v, u
        edges.add((u, v))

    G = {label: [] for label in labels}
    for u, v in edges:
        G[u].append(v)
        G[v].append(u)
    ##

    return G

if __name__ == '__main__':
    if len(sys.argv) != 7:
        raise SystemExit(f'Usage: {sys.argv[0]} <number of nodes> <number of edges> <minimum label length> <maximum label length> <alphabet size> <output_nodes_file> <output_edges_file>')
    n, m, min_label, max_label, alphabet_size, nodes_output, edges_output = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6], sys.argv[7]
    alphabet = [str(i) for i in range(alphabet_size)]

    G = generate_graph(n, m, min_label, max_label, alphabet)
    V_ids = {v: i+1 for i, v in enumerate(G.keys())}

    with open(nodes_output, 'w') as fout:
        print('left;right', file=fout)
        for i, node in enumerate(G.keys()):
            print(f'{V_ids[node]};{node}', file=fout)

    with open(edges_output, 'w') as fout:
        print('left;right', file=fout)
        for v in G:
            for u in G[v]:
                if v < u:
                    print(f'{V_ids[v]};{V_ids[u]}', file=fout)
