#!/usr/bin/env python3
"""
Generate a random graph instance with a target edge density.

Node labels are sequences sampled from a FASTA file. Edges are added
uniformly at random until the requested density is reached, seeded with an
initial random matching so that no vertex is isolated.

Usage:
    ./generate_instance.py --fasta sequences.fa --n 200 --density 0.05 \\
        --seed 42 --out-nodes nodes.fa --out-edges edges.int
"""

import argparse
import random
import sys
from fasta_utils import read_fasta

def add_edges(node_ids: list[str], density: float) -> dict[str, list[str]]:
    """
    Given a list of node identifiers, add random edges until `density` is
    reached. An initial random matching guarantees no isolated vertices.
    Edges are sampled uniformly at random.

    Returns an adjacency-list dict {node_id: [neighbour, ...]}.
    """
    num_nodes = len(node_ids)
    max_possible_edges = (num_nodes * (num_nodes - 1)) // 2
    num_edges = round(density * max_possible_edges)

    min_edges = (num_nodes + 1) // 2
    if num_edges < min_edges:
        raise ValueError(
            f'Requested density {density} yields {num_edges} edges, but at least '
            f'{min_edges} are needed to avoid isolated vertices (try a higher density).'
        )
    if num_edges > max_possible_edges:
        num_edges = max_possible_edges

    # Initial random matching; prevents isolated vertices.
    shuffled = node_ids[:]
    random.shuffle(shuffled)
    edges: set[tuple[str, str]] = set()
    for i in range(0, num_nodes - 1, 2):
        u, v = sorted([shuffled[i], shuffled[i + 1]])
        edges.add((u, v))
    if num_nodes % 2 == 1:  # odd node: connect to a random neighbour
        u, v = sorted([shuffled[-1], random.choice(shuffled[:-1])])
        edges.add((u, v))

    # Uniform random fill to reach target num_edges.
    while len(edges) < num_edges:
        u, v = sorted(random.sample(node_ids, 2))
        edges.add((u, v))

    adj_list: dict[str, list[str]] = {node_id: [] for node_id in node_ids}
    for u, v in sorted(edges):
        adj_list[u].append(v)
        adj_list[v].append(u)
    return adj_list

def make_fasta_graph(fasta_path: str, num_nodes: int,
                     density: float) -> tuple[dict[str, str], dict[str, list[str]]]:
    """Sample num_nodes sequences from a FASTA file then build a random graph."""
    sequences = read_fasta(fasta_path)
    if len(sequences) < num_nodes:
        raise ValueError(
            f'FASTA file contains only {len(sequences)} sequences but --n {num_nodes} '
            f'was requested.'
        )

    chosen_ids = random.sample(list(sequences.keys()), num_nodes)
    labels: dict[str, str] = {seq_id: sequences[seq_id] for seq_id in chosen_ids}

    adj_list = add_edges(chosen_ids, density)
    return labels, adj_list

def write_outputs(labels: dict[str, str],
                  adj_list: dict[str, list[str]],
                  out_nodes: str,
                  out_edges: str) -> None:
    # Assign stable numeric IDs in sorted order of the original node keys.
    node_ids = sorted(labels.keys())
    node_index = {node_id: i + 1 for i, node_id in enumerate(node_ids)}

    with open(out_nodes, 'w') as fout:
        for node_id in node_ids:
            fout.write(f'>{node_index[node_id]}\n')
            fout.write(f'{labels[node_id]}\n')

    with open(out_edges, 'w') as fout:
        for v in node_ids:
            for u in adj_list[v]:
                if v < u:
                    fout.write(f'{node_index[v]} {node_index[u]}\n')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('-f', '--fasta', required=True,
                        help='Input FASTA file to sample node labels from.')
    parser.add_argument('--n', type=int, required=True,
                        help='Number of nodes to sample.')
    parser.add_argument('--density', type=float, required=True,
                        help='Target edge density in [0, 1].')
    parser.add_argument('--seed', type=int, default=None,
                        help='Random seed for reproducibility.')
    parser.add_argument('-o', '--out-nodes', default='out_nodes.fa',
                        help='Output FASTA path.')
    parser.add_argument('-e', '--out-edges', default='out_edges.int',
                        help='Output interactions path.')

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
        labels, adj_list = make_fasta_graph(args.fasta, args.n, args.density)
        write_outputs(labels, adj_list, args.out_nodes, args.out_edges)
    except ValueError as exc:
        raise SystemExit(f'Error: {exc}') from None
    except OSError as exc:
        raise SystemExit(f'Error: {exc}') from None

    num_nodes = len(labels)
    num_edges = sum(len(neighbours) for neighbours in adj_list.values()) // 2
    density = num_edges / ((num_nodes * (num_nodes - 1)) // 2)
    sys.stderr.write(f'nodes={num_nodes}  edges={num_edges}  density={density:.4f}\n')

if __name__ == '__main__':
    main()
