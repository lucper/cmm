#!/usr/bin/env python3
"""Remove FASTA records whose sequence contains non-amino-acid letters,
drop interaction edges incident to any removed record, then remove nodes
left isolated (no incident edges).

Usage:
    ./clean_dataset.py --fasta sequences.fa --edges links.int \\
        --out-nodes sequences.clean.fa --out-edges links.clean.int
"""

import argparse
import sys
from fasta_utils import AA_ALPHABET_EXTENDED, parse_fasta

VALID_RESIDUES = set(AA_ALPHABET_EXTENDED)

def remove_dubious_nodes(fasta_in, edges_in, out_nodes, out_edges,
                         uppercase=False, keep_isolated=False):
    # Validate sequences. Buffer records that pass so we can later filter out
    # any that become isolated; remember insertion order.
    valid_records = []  # list of (header, seq_id, seq)
    valid_ids = set()
    num_seqs_total = 0

    for header, seq_id, seq in parse_fasta(fasta_in):
        num_seqs_total += 1
        if uppercase:
            seq = seq.upper()
        if seq and set(seq) <= VALID_RESIDUES:
            valid_records.append((header, seq_id, seq))
            valid_ids.add(seq_id)

    # Keep edges whose both endpoints have a valid sequence, and record which
    # nodes still have at least one incident edge.
    num_edges_total = num_edges_kept = 0
    non_isolated = set()
    with open(edges_in) as fh, open(out_edges, "w") as fout:
        for line in fh:
            if not line.strip():
                continue
            num_edges_total += 1
            fields = line.split()
            if len(fields) < 2:
                sys.stderr.write(f"Warning: skipping malformed edge line: {line!r}\n")
                continue
            u, v = fields[0], fields[1]
            if u in valid_ids and v in valid_ids:
                fout.write(f"{u} {v}\n")
                num_edges_kept += 1
                non_isolated.add(u)
                non_isolated.add(v)

    # Write FASTA records, dropping isolated ones unless asked to keep.
    num_seqs_kept = 0
    with open(out_nodes, "w") as fout:
        for header, seq_id, seq in valid_records:
            if keep_isolated or seq_id in non_isolated:
                fout.write(f"{header}\n{seq}\n")
                num_seqs_kept += 1

    num_invalid_removed = num_seqs_total - len(valid_records)
    num_isolated_removed = len(valid_records) - num_seqs_kept
    sys.stderr.write(
        f"Kept sequences:  {num_seqs_kept}  (was {num_seqs_total}; "
        f"removed {num_invalid_removed} invalid, {num_isolated_removed} isolated)\n"
        f"Kept edges:      {num_edges_kept}  (was {num_edges_total})\n"
    )

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('-f', '--fasta', required=True,
                        help='Input FASTA file with protein sequences.')
    parser.add_argument('-i', '--edges', required=True,
                        help='Input interactions file with lines formatted as "u v".')
    parser.add_argument('-o', '--out-nodes', default=None,
                        help='Output FASTA path (default: <fasta>.clean).')
    parser.add_argument('-e', '--out-edges', default=None,
                        help='Output interactions path (default: <edges>.clean).')
    parser.add_argument('-u', '--uppercase', action='store_true',
                        help='Uppercase sequences before validating (keeps soft-masked residues).')
    parser.add_argument('--keep-isolated', action='store_true',
                        help='Keep valid sequences even if they have no incident edges '
                             '(default: remove isolated nodes).')

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    out_nodes = args.out_nodes or args.fasta + '.clean'
    out_edges = args.out_edges or args.edges + '.clean'

    try:
        remove_dubious_nodes(args.fasta, args.edges, out_nodes, out_edges,
                             args.uppercase, args.keep_isolated)
    except OSError as exc:
        raise SystemExit(f'Error: {exc}') from None

if __name__ == '__main__':
    main()
