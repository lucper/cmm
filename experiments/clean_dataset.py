#!/usr/bin/env python3
"""Remove FASTA records whose sequence contains non-amino-acid letters,
drop interaction edges incident to any removed record, then remove nodes
left isolated (no incident edges).
"""

import argparse
import sys

AA_ALPHABET = set("ARNDCEQGHILKMFPSTWYVU")

def parse_fasta(path):
    """Yield (header_line, seq_id, sequence) for each record.

    seq_id matches how the importer reads it: text after '>', truncated at
    the first whitespace. Sequence lines are concatenated with whitespace
    stripped, mirroring how wrapped FASTA lines are joined.
    """
    header = seq_id = None
    chunks = []

    def emit():
        if header is not None:
            yield header, seq_id, "".join(chunks)

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                yield from emit()
                header = line
                seq_id = line[1:].split(None, 1)[0] if line[1:].split() else ""
                chunks = []
            elif line.strip():
                chunks.append("".join(line.split()))
    yield from emit()

def remove_dubious_nodes(fasta_in, int_in, fasta_out, int_out, uppercase=False,
          keep_isolated=False):
    # Validate sequences. Buffer records that pass so we can later filter out
    # any that become isolated; remember insertion order.
    valid_records = [] # list of (header, seq_id, seq)
    valid_ids = set()
    seqs_total = 0

    for header, seq_id, seq in parse_fasta(fasta_in):
        seqs_total += 1
        if uppercase:
            seq = seq.upper()
        if seq and set(seq) <= AA_ALPHABET:
            valid_records.append((header, seq_id, seq))
            valid_ids.add(seq_id)

    # Keep edges whose both endpoints have a valid sequence, and record which
    # nodes still have at least one incident edge.
    edges_total = edges_kept = 0
    non_isolated = set()
    with open(int_in) as fh, open(int_out, "w") as out:
        for line in fh:
            if not line.strip():
                continue
            edges_total += 1
            parts = line.split()
            if len(parts) < 2:
                sys.stderr.write(f"Warning: skipping malformed edge line: {line!r}\n")
                continue
            u, v = parts[0], parts[1]
            if u in valid_ids and v in valid_ids:
                out.write(f"{u} {v}\n")
                edges_kept += 1
                non_isolated.add(u)
                non_isolated.add(v)

    # Write FASTA records, dropping isolated ones unless asked to keep.
    seqs_kept = 0
    with open(fasta_out, "w") as out:
        for header, seq_id, seq in valid_records:
            if keep_isolated or seq_id in non_isolated:
                out.write(f"{header}\n{seq}\n")
                seqs_kept += 1

    invalid_removed = seqs_total - len(valid_records)
    isolated_removed = len(valid_records) - seqs_kept
    sys.stderr.write(
        f"Kept sequences:  {seqs_kept}  (was {seqs_total}; "
        f"removed {invalid_removed} invalid, {isolated_removed} isolated)\n"
        f"Kept edges:      {edges_kept}  (was {edges_total})\n"
    )

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fasta_in")
    p.add_argument("int_in")
    p.add_argument("-o", "--fasta-out", default=None,
                   help="cleaned FASTA output (default: <fasta_in>.clean)")
    p.add_argument("-e", "--int-out", default=None,
                   help="cleaned interaction output (default: <int_in>.clean)")
    p.add_argument("-u", "--uppercase", action="store_true",
                   help="uppercase sequences before validating (keeps soft-masked residues)")
    p.add_argument("--keep-isolated", action="store_true",
                   help="keep valid sequences even if they have no incident edges "
                        "(default: remove isolated nodes)")
    args = p.parse_args()

    fasta_out = args.fasta_out or args.fasta_in + ".clean"
    int_out = args.int_out or args.int_in + ".clean"
    remove_dubious_nodes(args.fasta_in, args.int_in, fasta_out, int_out, args.uppercase,
          args.keep_isolated)

if __name__ == "__main__":
    main()
