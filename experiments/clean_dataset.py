#!/usr/bin/env python3
"""Remove FASTA records whose sequence contains non-amino-acid letters,
then drop interaction edges incident to any removed record.
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

def clean(fasta_in, int_in, fasta_out, int_out, uppercase=False):
    kept = set()
    seqs_total = 0

    with open(fasta_out, "w") as out:
        for header, seq_id, seq in parse_fasta(fasta_in):
            seqs_total += 1
            if uppercase:
                seq = seq.upper()
            if seq and set(seq) <= AA_ALPHABET:
                kept.add(seq_id)
                out.write(f"{header}\n{seq}\n")

    edges_total = edges_kept = 0
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
            if u in kept and v in kept:
                out.write(f"{u} {v}\n")
                edges_kept += 1

    sys.stderr.write(
        f"Kept sequences:  {len(kept)}  (was {seqs_total})\n"
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
    args = p.parse_args()

    fasta_out = args.fasta_out or args.fasta_in + ".clean"
    int_out = args.int_out or args.int_in + ".clean"
    clean(args.fasta_in, args.int_in, fasta_out, int_out, args.uppercase)

if __name__ == "__main__":
    main()
