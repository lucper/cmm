#!/usr/bin/env python3

import argparse
import sys

WILDCARD = 'x'

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

def read_motif_pairs(path: str) -> list[tuple[str, str]]:
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
            if len(parts) < 3:
                print(
                    f"  [WARN] truth file line {lineno}: fewer than 3 columns — skipped.",
                    file=sys.stderr,
                )
                continue
            pairs.append((parts[0], parts[1]) if parts[0] <= parts[1] else (parts[1], parts[0]))
    return pairs

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

def similarity(pair_A: tuple[str, str], pair_B: tuple[str, str], V: dict[str, str], E: list[tuple[str, str]]) -> float:
    """
    Returns
        s(pair_A, pair_B, G=(V,E)) = pos(E_A, E_B) / |E_A \\cup E_B|,
    where E_A (resp. E_B) is the number of edges in G in which the motif pair `pair_A`
    (resp. `pair_B`) co-occur and pos(E_A, E_B) is the number of edges {u,v} in which
    there exists substrings u' of u and v' of v such that
        - u' matches pair_A[0] and pair_B[0]; and
        - v' matches pair_A[1] and pair_B[1].
    """
    pos = {motif : {nid: positions
                    for nid, seq in V.items()
                    if (positions := find_motif_positions(motif, seq))}
           for motif in {*pair_A, *pair_B}}

    E_A, E_B, pos_AB = set(), set(), set()

    for u, v in E:
        u, v = (u, v) if u <= v else (v, u)

        a0_u = pos[pair_A[0]].get(u); a1_v = pos[pair_A[1]].get(v)
        a0_v = pos[pair_A[0]].get(v); a1_u = pos[pair_A[1]].get(u)
        pair_A_occurs = (a0_u and a1_v) or (a0_v and a1_u)

        b0_u = pos[pair_B[0]].get(u); b1_v = pos[pair_B[1]].get(v)
        b0_v = pos[pair_B[0]].get(v); b1_u = pos[pair_B[1]].get(u)
        pair_B_occurs = (b0_u and b1_v) or (b0_v and b1_u)

        if pair_A_occurs: E_A.add((u, v))
        if pair_B_occurs: E_B.add((u, v))

        if ((a0_u and b0_u and (a0_u & b0_u)) and
            (a1_v and b1_v and (a1_v & b1_v))):
            pos_AB.add((u, v))
        elif ((a0_v and b0_v and (a0_v & b0_v)) and
              (a1_u and b1_u and (a1_u & b1_u))):
            pos_AB.add((u, v))

    return len(pos_AB) / len(E_A | E_B) if E_A or E_B else 0.0

def sensitivity(S: list[tuple[str, str]], T: list[tuple[str, str]], V: dict[str, str], E: list[tuple[str, str]], k: int) -> float:
    if not S or not T:
        return 0.0

    S_star = []
    for candidate in S:
        if not any(similarity(candidate, kept_pair, V, E) == 1.0
                   for kept_pair in S_star):
            S_star.append(candidate)

    if k < 1 or k > len(S_star):
        raise ValueError(f"k={k} must satisfy 1 <= k <= |S*|={len(S_star)}")

    tp_count = sum(
        1
        for pred in S_star[:k]
        if any(
            similarity(pred, truth_pair, V, E) == 1.0
            for truth_pair in T
        )
    )

    return tp_count / k

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
        help="Output file produced by motif mining algorithm "
             "(3 columns: motif_X motif_Y f_score)"
    )
    parser.add_argument(
        "--k", required=True, type=int, default=1,
        help="Number of 'predicted' motif pairs to compare against 'truth'"
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

    print(f"[1/4] Reading truth set     : {args.truth}", file=sys.stderr)
    T = read_motif_pairs(args.truth)
    if not T:
        raise SystemExit(f"{args.truth} is empty")
    print(f"        {len(T)} truth pair(s) loaded.", file=sys.stderr)

    print(f"[2/4] Reading predictions   : {args.sequences}", file=sys.stderr)
    S = read_motif_pairs(args.predicted)
    if not S:
        raise SystemExit(f"{args.predicted} is empty.")
    print(f"        {len(S)} sequences loaded.", file=sys.stderr)

    print(f"[3/4] Reading sequences     : {args.sequences}", file=sys.stderr)
    V = read_fasta(args.sequences)
    print(f"        {len(V)} sequence(s) loaded.", file=sys.stderr)

    print(f"[4/4] Reading interactions  : {args.interactions}", file=sys.stderr)
    E = read_interactions(args.interactions)
    print(f"        {len(E)} interaction(s) loaded.", file=sys.stderr)

    print(f"\nComputing sensitivity for k={args.k}...\n", file=sys.stderr)
    result = sensitivity(S, T, V, E, args.k)
    print(f"|T| (truth pairs)       : {len(T):>6}")
    print(f"|S| (predicted pairs)   : {len(S):>6}")
    print(f" k  (cutoff)            : {args.k:>6}")
    print(f"{'-'*35}")
    print(f"Sensitivity(S, T, k)    : {result:.3f}")

if __name__ == "__main__":
    main()
