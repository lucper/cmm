#!/usr/bin/env python3

# The 20 standard amino acids, sorted. Used for synthetic label generation,
# where the alphabet is sliced by size.
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

# Standard residues plus selenocysteine (U). Used when validating sequences
# coming from real databases such as STRING.
AA_ALPHABET_EXTENDED = AA_ALPHABET + "U"

def parse_fasta(path):
    """Yield (header_line, seq_id, sequence) for each record.

    seq_id is the text after '>' truncated at the first whitespace, matching
    how read_graph_files() in data_import.cpp reads it. Sequence lines are
    concatenated with all whitespace stripped, mirroring how wrapped FASTA
    lines are joined by the importer.
    """
    header = seq_id = None
    chunks = []

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, seq_id, "".join(chunks)
                header = line
                fields = line[1:].split(None, 1)
                seq_id = fields[0] if fields else ""
                chunks = []
            elif line.strip():
                chunks.append("".join(line.split()))

    if header is not None:
        yield header, seq_id, "".join(chunks)


def read_fasta(path):
    """Return {seq_id: sequence} preserving insertion order, uppercased."""
    return {seq_id: seq.upper() for _, seq_id, seq in parse_fasta(path)}
