#!/usr/bin/python3
"""
Generate a random graph instance with a target edge density.

Two modes:
  Synthetic  (default): labels are randomly generated strings.
  FASTA      (--fasta): labels are sequences sampled from a FASTA file.

In both modes edges are added uniformly at random until the requested
density is reached, seeded with an initial random matching so that no
vertex is isolated.

Usage
synthetic mode:
    python3 generate_instance_density.py \
        --n 200 --density 0.05 \
        --min-label-len 4 --max-label-len 8 --alphabet-size 20 \
        --seed 42 \
        --out-nodes nodes.fasta --out-edges edges.txt

FASTA mode:
    python3 generate_instance_density.py \
        --fasta sequences.fasta --n 200 --density 0.05 \
        --seed 42 \
        --out-nodes nodes.fasta --out-edges edges.txt
"""

import argparse
import random
import string
import sys

AA_ALPHABET = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def read_fasta(path: str) -> dict[str, str]:
    """Return {id: sequence} preserving insertion order."""
    sequences: dict[str, str] = {}
    current_id = None
    buf: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(buf)
                current_id = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.upper())
    if current_id is not None:
        sequences[current_id] = ''.join(buf)
    return sequences

def add_edges(node_ids: list[str], density: float) -> dict[str, list[str]]:
    """
    Given a list of node identifiers, add random edges until `density` is
    reached.  An initial random matching guarantees no isolated vertices.
    Edges are sampled uniformly at random.

    Returns an adjacency-list dict {node_id: [neighbour, ...]}.
    """
    n = len(node_ids)
    max_possible_edges = (n * (n - 1)) // 2
    m = round(density * max_possible_edges)

    min_m = (n + 1) // 2
    if m < min_m:
        raise ValueError(
            f'Requested density {density} yields {m} edges, but at least {min_m} '
            f'are needed to avoid isolated vertices (try a higher density).'
        )
    if m > max_possible_edges:
        m = max_possible_edges

    # Initial random matching - prevents isolated vertices without imposing structure.
    shuffled = node_ids[:]
    random.shuffle(shuffled)
    edges: set[tuple[str, str]] = set()
    for i in range(0, n - 1, 2):
        u, v = sorted([shuffled[i], shuffled[i + 1]])
        edges.add((u, v))
    if n % 2 == 1:                          # odd node: connect to a random neighbour
        u, v = sorted([shuffled[-1], random.choice(shuffled[:-1])])
        edges.add((u, v))

    # Uniform random fill to reach target m.
    while len(edges) < m:
        u, v = sorted(random.sample(node_ids, 2))
        edges.add((u, v))

    G: dict[str, list[str]] = {nid: [] for nid in node_ids}
    for u, v in sorted(edges):
        G[u].append(v)
        G[v].append(u)
    return G

def make_synthetic_graph(n: int, density: float,
                         min_label_len: int, max_label_len: int,
                         alphabet: str) -> tuple[dict[str, str], dict[str, list[str]]]:
    """Generate n random labels then build a random graph."""
    max_possible_labels = sum(
        len(alphabet) ** l for l in range(min_label_len, max_label_len + 1)
    )
    if n > max_possible_labels:
        raise ValueError(
            f'Too many nodes: max is {max_possible_labels} for the given alphabet and label lengths.'
        )

    label_set: set[str] = set()
    while len(label_set) < n:
        length = random.randint(min_label_len, max_label_len)
        label_set.add(''.join(random.choices(alphabet, k=length)))

    node_ids = sorted(label_set)
    labels: dict[str, str] = {nid: nid for nid in node_ids}

    G = add_edges(node_ids, density)
    return labels, G

def make_fasta_graph(fasta_path: str, n: int,
                     density: float) -> tuple[dict[str, str], dict[str, list[str]]]:
    """Sample n sequences from a FASTA file then build a random graph."""
    all_seqs = read_fasta(fasta_path)
    if len(all_seqs) < n:
        raise ValueError(
            f'FASTA file contains only {len(all_seqs)} sequences but --n {n} was requested.'
        )

    chosen_ids = random.sample(list(all_seqs.keys()), n)
    labels: dict[str, str] = {sid: all_seqs[sid] for sid in chosen_ids}

    G = add_edges(chosen_ids, density)
    return labels, G

def write_outputs(labels: dict[str, str],
                  G: dict[str, list[str]],
                  nodes_output: str,
                  edges_output: str) -> None:
    # Assign stable numeric IDs in sorted order of the original node keys.
    node_ids = sorted(labels.keys())
    V_ids = {nid: i + 1 for i, nid in enumerate(node_ids)}

    with open(nodes_output, 'w') as fout:
        for nid in node_ids:
            print(f'>{V_ids[nid]}', file=fout)
            print(labels[nid], file=fout)

    with open(edges_output, 'w') as fout:
        for v in node_ids:
            for u in G[v]:
                if v < u:
                    print(f'{V_ids[v]} {V_ids[u]}', file=fout)

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            'Generate a random graph instance with a target edge density.\n\n'
            'Two modes:\n'
            '  Synthetic (default)  - labels are randomly generated strings.\n'
            '  FASTA    (--fasta)   - labels are sequences sampled from a FASTA file.\n\n'
            'In both modes edges are added uniformly at random until the\n'
            'requested density is reached.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('--n',        type=int,   required=True,  help='Number of nodes to sample/generate.')
    parser.add_argument('--density',  type=float, required=True,  help='Target edge density in [0, 1].')
    parser.add_argument('--seed',     type=int,   default=None,   help='Random seed for reproducibility.')
    parser.add_argument('--out-nodes', default='out_nodes.fasta', help='Output FASTA path.')
    parser.add_argument('--out-edges', default='out_edges.txt',   help='Output interactions path.')

    # FASTA mode
    fasta_group = parser.add_argument_group('FASTA mode (mutually exclusive with synthetic mode)')
    fasta_group.add_argument('--fasta', default=None,
                             help='Input FASTA file; when given, node labels are sampled from it.')

    # Synthetic mode
    syn_group = parser.add_argument_group('Synthetic mode (ignored when --fasta is given)')
    syn_group.add_argument('--min-label-len', type=int, default=4,  help='Minimum label length.')
    syn_group.add_argument('--max-label-len', type=int, default=8,  help='Maximum label length.')
    syn_group.add_argument('--alphabet-size', type=int, default=20, help='Alphabet size (up to 20).')

    return parser.parse_args()

def main() -> None:
    args = parse_args()

    if not (0.0 <= args.density <= 1.0):
        raise SystemExit('--density must be in [0.0, 1.0].')
    if args.n < 2:
        raise SystemExit('--n must be at least 2.')

    if args.seed is not None:
        random.seed(args.seed)

    try:
        if args.fasta:
            # FASTA mode
            labels, G = make_fasta_graph(args.fasta, args.n, args.density)
        else:
            # Synthetic mode
            ascii_chars = string.ascii_lowercase + string.ascii_uppercase + string.digits
            if args.alphabet_size > len(ascii_chars):
                raise SystemExit(
                    f'--alphabet-size ({args.alphabet_size}) must be at most {len(AA_ALPHABET)}.'
                )
            if args.min_label_len > args.max_label_len:
                raise SystemExit('--min-label-len must be <= --max-label-len.')
            alphabet = AA_ALPHABET[:args.alphabet_size]
            labels, G = make_synthetic_graph(
                args.n, args.density,
                args.min_label_len, args.max_label_len,
                alphabet,
            )
    except ValueError as exc:
        raise SystemExit(f'Error: {exc}') from None

    write_outputs(labels, G, args.out_nodes, args.out_edges)

    n = len(labels)
    actual_edges = sum(len(v) for v in G.values()) // 2
    actual_density = actual_edges / ((n * (n - 1)) // 2)
    print(f'nodes={n}  edges={actual_edges}  density={actual_density:.4f}', file=sys.stderr)


if __name__ == '__main__':
    main()
