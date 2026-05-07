#!/usr/bin/env python3

import argparse
import sys
from pprint import pprint

WILDCARD = 'x'

def standardise(a: str, b: str) -> tuple[str, str]:
    """Return (min, max) under lexicographic order so X <= Y always holds."""
    return (a, b) if a <= b else (b, a)

def read_fasta(path: str) -> dict[str, str]:
    """Return OrderedDict {id: sequence} from a FASTA file."""
    sequences: dict[str, str] = {}
    current_id = None
    buf: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(buf)
                # Use the first word (before a space) in the header as ID
                current_id = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.upper())
    if current_id is not None:
        sequences[current_id] = "".join(buf)
    return sequences

def read_interactions(path: str) -> list[tuple[str, str]]:
    """Return list of (id1, id2) pairs from a space-separated file."""
    edges: list[tuple[str, str]] = []
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) == 2:
                edges.append((parts[0], parts[1]))
            else:
                raise ValueError(f"Illegal line in {path}: {line.strip()}")
    return edges

def read_truth(path: str) -> list[tuple[str, str]]:
    """
    Load the truth-set file (3 columns: motif_X motif_Y score).
    Returns a set of standardised (X, Y) pairs.
    """
    pairs: list[tuple[str, str]] = list()
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                print(
                    f"  [WARN] truth file line {lineno}: fewer than 2 columns — skipped.",
                    file=sys.stderr,
                )
                continue
            pairs.append(standardise(parts[0], parts[1]))
    return pairs

def read_predictions(path: str) -> list[tuple[str, str]]:
    """
    Load the mining algorithm output (3 columns: motif_A motif_B score).
    Returns a list of (motif_A, motif_B, score_str) with motifs upper-cased.
    Lines with fewer than 3 columns are warned and skipped.
    """
    rows: list[tuple[str, str, str]] = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                print(
                    f"  [WARN] prediction file line {lineno}: fewer than 3 columns — skipped.",
                    file=sys.stderr,
                )
                continue
            rows.append((parts[0], parts[1]))
    return rows

def motif_matches(motif: str, text: str, start: int) -> bool:
    """Check whether `motif` matches `text` starting at `start`."""
    for i, ch in enumerate(motif):
        if ch != WILDCARD and text[start + i] != ch:
            return False
    return True

def find_motif_positions(motif: str, sequence: str) -> set[int]:
    """Return all start positions where `motif` matches `sequence`."""
    l = len(motif)
    return {i for i in range(len(sequence) - l + 1) if motif_matches(motif, sequence, i)}

def similarity(pair_A: tuple[str, str], pair_B: tuple[str, str], V: dict[str, str], E: list[tuple[str, str]]):
    """
    Returns
        s(pair_A, pair_B, G=(V,E)) = pos(E_A, E_B) / |E_A \\cup E_B|,
    where E_A (resp. E_B) is the number of edges in G in which the motif pair `pair_A`
    (resp. `pair_B`) co-occur and pos(E_A, E_B) is the number of edges {u,v} in which
    there exists substrings u' of u and v' of v such that
        - u' matches pair_A[0] and pair_B[0]; and
        - v' matches pair_A[1] and pair_B[1].
    """
    E_ordered = [standardise(u, v) for u, v in E] # necessary?
    E_A = {(u, v) for u, v in E_ordered if find_motif_positions(pair_A[0], V[u]) and find_motif_positions(pair_A[1], V[v])}
    E_B = {(u, v) for u, v in E_ordered if find_motif_positions(pair_B[0], V[u]) and find_motif_positions(pair_B[1], V[v])}
    pos_AB = {(u, v) for u, v in E_ordered if (find_motif_positions(pair_A[0], V[u]) & find_motif_positions(pair_B[0], V[u])) and
                                              (find_motif_positions(pair_A[1], V[v]) & find_motif_positions(pair_B[1], V[v]))}
    return len(pos_AB) / len(E_A | E_B)

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate motif-pair mining output against a truth set. "
            "Pairs are standardised to (X,Y) with X<=Y (lex) before comparison."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--truth", required=True, default=argparse.SUPPRESS,
        help="Truth-set file produced by ppi_network_generator.py "
             "(3 columns: motif_X motif_Y f_score)"
    )
    parser.add_argument(
        "--predicted", required=True, default=argparse.SUPPRESS,
        help="Mining algorithm output (3 columns: motif_A motif_B score)"
    )
    parser.add_argument(
        "--k", required=True, type=int, default=1,
        help="Rank to compare 'truth' against 'predicted'"
    )
    parser.add_argument(
        "--sequences", required=True, default=argparse.SUPPRESS,
        help="FASTA file with sequences"
    )
    parser.add_argument(
        "--interactions", required=True, default=argparse.SUPPRESS,
        help="Text file with two space-separated columns specifying interactions between sequences"
    )
    return parser.parse_args()

def main() -> None:
    args = parse_args()

    truth = read_truth(args.truth)
    predictions = read_predictions(args.predicted)
    nodes = read_fasta(args.sequences)
    edges = read_interactions(args.interactions)
    
    # ----- Example -----
    pair_A, pair_B = truth[0], predictions[4]
    similarity(pair_A, pair_B, nodes, edges)

if __name__ == "__main__":
    main()
