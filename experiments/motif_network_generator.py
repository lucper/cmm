#!/usr/bin/env python3
"""
Protein Interaction Network Motif Generator
============================================
Generates a synthetic protein interaction dataset with implanted (l,d)-motifs.

Usage:
    python motif_network_generator.py \\
        --fasta proteins.fasta \\
        --interactions interactions.txt \\
        --l 10 --d 2 --K 5 --S 50 --E 0.3 \\
        [--seed 42] \\
        [--out-fasta out.fasta] \\
        [--out-interactions out_interactions.txt] \\
        [--out-motifs out_motifs.txt]
"""

import argparse
import random
import sys
from pathlib import Path
from itertools import combinations


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

AA_ALPHABET = list("ACDEFGHIKLMNPQRSTVWY")
WILDCARD = "."          # regex-style wildcard stored in motif strings
NUM_CANDIDATE_MOTIFS = 50


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

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
                current_id = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.upper())
    if current_id is not None:
        sequences[current_id] = "".join(buf)
    return sequences


def write_fasta(path: str, sequences: dict[str, str]) -> None:
    with open(path, "w") as fh:
        for sid, seq in sequences.items():
            fh.write(f">{sid}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")


def read_interactions(path: str) -> list[tuple[str, str]]:
    """Return list of (id1, id2) pairs from a space-separated file."""
    edges: list[tuple[str, str]] = []
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) >= 2:
                edges.append((parts[0], parts[1]))
    return edges


def write_interactions(path: str, edges: set[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for u, v in sorted(edges):
            fh.write(f"{u} {v}\n")


def write_motifs(path: str, motifs: list[str]) -> None:
    with open(path, "w") as fh:
        for i, m in enumerate(motifs, 1):
            fh.write(f"motif_{i:03d}\t{m}\n")


# ---------------------------------------------------------------------------
# (l,d)-motif helpers
# ---------------------------------------------------------------------------

def generate_random_motif(l: int, d: int) -> str:
    """
    Generate one random (l,d)-motif: length l, exactly d wildcard positions.
    Non-wildcard positions are sampled from AA_ALPHABET.
    """
    if d > l:
        raise ValueError(f"d={d} cannot exceed l={l}")
    positions = list(range(l))
    wild_positions = set(random.sample(positions, d))
    return "".join(
        WILDCARD if i in wild_positions else random.choice(AA_ALPHABET)
        for i in range(l)
    )


def motif_matches(motif: str, text: str, start: int) -> bool:
    """Check whether `motif` matches `text` starting at `start`."""
    for i, ch in enumerate(motif):
        if ch != WILDCARD and text[start + i] != ch:
            return False
    return True


def find_motif_positions(motif: str, sequence: str) -> list[int]:
    """Return all start positions where `motif` matches `sequence`."""
    l = len(motif)
    return [i for i in range(len(sequence) - l + 1) if motif_matches(motif, sequence, i)]


def sequence_contains_motif(motif: str, sequence: str) -> bool:
    return bool(find_motif_positions(motif, sequence))


def implant_motif(sequence: str, motif: str) -> str:
    """
    Implant `motif` at a random valid position in `sequence`.
    Wildcard positions keep the original amino acid; fixed positions are overwritten.
    Raises ValueError if the sequence is shorter than the motif.
    """
    l = len(motif)
    if len(sequence) < l:
        raise ValueError("Sequence shorter than motif length.")
    start = random.randint(0, len(sequence) - l)
    seq_list = list(sequence)
    for i, ch in enumerate(motif):
        if ch != WILDCARD:
            seq_list[start + i] = ch
    return "".join(seq_list)


# ---------------------------------------------------------------------------
# Edge-density helper
# ---------------------------------------------------------------------------

def edge_density(edges: set[tuple[str, str]], n: int) -> float:
    """Fraction of possible undirected edges that are present (max = n*(n-1)/2)."""
    max_edges = n * (n - 1) / 2
    if max_edges == 0:
        return 0.0
    return len(edges) / max_edges


# ---------------------------------------------------------------------------
# Core pipeline
# ---------------------------------------------------------------------------

def run(
    fasta_path: str,
    interactions_path: str,
    l: int,
    d: int,
    K: int,
    S: int,
    E: float,
    seed: int | None,
    out_fasta: str,
    out_interactions: str,
    out_motifs: str,
) -> None:

    # -- Reproducibility ---------------------------------------------------
    if seed is not None:
        random.seed(seed)

    # -- Load inputs -------------------------------------------------------
    print("[1/7] Reading inputs …")
    all_seqs = read_fasta(fasta_path)
    all_edges_raw = read_interactions(interactions_path)
    all_ids = list(all_seqs.keys())

    if len(all_ids) < S:
        sys.exit(
            f"ERROR: FASTA contains only {len(all_ids)} sequences but S={S} was requested."
        )

    # -- Sample S sequences ------------------------------------------------
    print(f"[2/7] Sampling S={S} sequences …")
    chosen_ids = random.sample(all_ids, S)
    chosen_set = set(chosen_ids)
    sequences: dict[str, str] = {sid: all_seqs[sid] for sid in chosen_ids}

    # Subset interactions to chosen proteins
    original_edges: set[tuple[str, str]] = set()
    for u, v in all_edges_raw:
        if u in chosen_set and v in chosen_set:
            a, b = (u, v) if u <= v else (v, u)
            original_edges.add((a, b))

    # -- Generate 50 random (l,d)-motifs, keep K of them ------------------
    print(f"[3/7] Generating {NUM_CANDIDATE_MOTIFS} candidate (l={l},d={d})-motifs …")
    candidate_motifs = [generate_random_motif(l, d) for _ in range(NUM_CANDIDATE_MOTIFS)]
    # Deduplicate (unlikely collision but defensive)
    candidate_motifs = list(dict.fromkeys(candidate_motifs))

    # -- Implant instances of ALL candidate motifs into sequences ----------
    print(f"[4/7] Implanting motif instances (I ∈ [3,10] per motif) …")
    # motif_hosts[motif_idx] = set of protein IDs that contain that motif
    motif_hosts: dict[int, set[str]] = {i: set() for i in range(len(candidate_motifs))}

    for motif_idx, motif in enumerate(candidate_motifs):
        I = random.randint(3, 10)               # instances for this motif
        # Shuffle IDs and pick I unique proteins
        targets = random.sample(chosen_ids, min(I, S))
        for sid in targets:
            sequences[sid] = implant_motif(sequences[sid], motif)
            motif_hosts[motif_idx].add(sid)

    # -- Build a "perfect" network by pairing motifs until density ≥ E ----
    print(f"[5/7] Building perfect interaction network (target density E={E}) …")
    implanted_edges: set[tuple[str, str]] = set()

    # All possible motif pairs
    motif_indices = list(range(len(candidate_motifs)))
    all_motif_pairs = list(combinations(motif_indices, 2))
    random.shuffle(all_motif_pairs)

    used_motif_pairs: list[tuple[int, int]] = []

    for mx, my in all_motif_pairs:
        hosts_x = motif_hosts[mx]
        hosts_y = motif_hosts[my]
        # Add edges between all proteins containing X and all containing Y
        for u in hosts_x:
            for v in hosts_y:
                if u != v:
                    a, b = (u, v) if u <= v else (v, u)
                    implanted_edges.add((a, b))
        used_motif_pairs.append((mx, my))
        current_density = edge_density(implanted_edges, S)
        if current_density >= E:
            print(
                f"    Density {current_density:.4f} reached after "
                f"{len(used_motif_pairs)} motif pair(s)."
            )
            break
    else:
        current_density = edge_density(implanted_edges, S)
        print(
            f"    WARNING: All motif pairs exhausted. "
            f"Final density = {current_density:.4f} (target was {E})."
        )

    # -- "Perfect" network: interactions ↔ motif-pair presence -------------
    # The implanted_edges IS the perfect network.  We discard all original
    # edges and replace them entirely so the result is perfectly explained
    # by the implanted motif pairs.
    final_edges = implanted_edges

    # -- Verify motifs actually present ------------------------------------
    # After ALL implantations are done, confirm each motif is still detectable.
    # (Earlier implants can be overwritten by later ones at overlapping sites,
    # so we re-implant any that went missing before verifying.)
    print("[6/7] Verifying / re-implanting any overwritten motifs …")
    repaired = 0
    active_motif_indices = list(set(i for pair in used_motif_pairs for i in pair))

    # Repair pass: iterate motif-by-motif and protein-by-protein.
    # To avoid circular overwrite (motif A repairs then motif B overwrites A),
    # we collect ALL (motif, protein) slots that are missing and implant them
    # in a single sequential sweep, checking after each implant.
    missing_slots = [
        (motif_idx, sid)
        for motif_idx in active_motif_indices
        for sid in motif_hosts[motif_idx]
        if not sequence_contains_motif(candidate_motifs[motif_idx], sequences[sid])
    ]

    MAX_REPAIR_ROUNDS = 10
    for _ in range(MAX_REPAIR_ROUNDS):
        if not missing_slots:
            break
        for motif_idx, sid in missing_slots:
            sequences[sid] = implant_motif(sequences[sid], candidate_motifs[motif_idx])
            repaired += 1
        # Recompute missing after this round
        missing_slots = [
            (motif_idx, sid)
            for motif_idx in active_motif_indices
            for sid in motif_hosts[motif_idx]
            if not sequence_contains_motif(candidate_motifs[motif_idx], sequences[sid])
        ]

    # Final check
    mismatches = sum(
        1
        for motif_idx in active_motif_indices
        for sid in motif_hosts[motif_idx]
        if not sequence_contains_motif(candidate_motifs[motif_idx], sequences[sid])
    )

    if repaired:
        print(f"    Re-implanted {repaired} overwritten instance(s).")
    if mismatches:
        print(f"    WARNING: {mismatches} motif instance(s) still missing after repair!")
    else:
        print("    All motif instances verified ✓")

    # -- Select K output motifs (from the used ones) -----------------------
    unique_used_indices = list(dict.fromkeys(i for pair in used_motif_pairs for i in pair))
    if len(unique_used_indices) < K:
        print(
            f"    NOTE: Only {len(unique_used_indices)} distinct motifs were used to "
            f"reach density E; returning all of them (K={K} requested)."
        )
        output_motif_indices = unique_used_indices
    else:
        output_motif_indices = random.sample(unique_used_indices, K)
    output_motifs = [candidate_motifs[i] for i in output_motif_indices]

    # -- Write outputs -----------------------------------------------------
    print("[7/7] Writing outputs …")
    write_fasta(out_fasta, sequences)
    write_interactions(out_interactions, final_edges)
    write_motifs(out_motifs, output_motifs)

    # -- Summary -----------------------------------------------------------
    print("\n=== Summary ===")
    print(f"  Output FASTA          : {out_fasta}  ({S} sequences)")
    print(f"  Output interactions   : {out_interactions}  ({len(final_edges)} edges)")
    print(f"  Output motifs         : {out_motifs}  ({len(output_motifs)} motifs)")
    print(f"  Final edge density    : {edge_density(final_edges, S):.4f}")
    print(f"  Motif pairs used      : {len(used_motif_pairs)}")
    print(f"  (l, d)                : ({l}, {d})")
    print(f"  K (output motifs)     : {K}")
    print(f"  S (proteins)          : {S}")
    print(f"  E (density target)    : {E}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a synthetic protein interaction dataset with (l,d)-motifs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--fasta",         required=True,  help="Input FASTA file of protein sequences")
    parser.add_argument("--interactions",  required=True,  help="Input space-separated interaction pairs")
    parser.add_argument("--l",   type=int, required=True,  help="Motif length")
    parser.add_argument("--d",   type=int, required=True,  help="Number of wildcard positions per motif")
    parser.add_argument("--K",   type=int, required=True,  help="Number of motifs to include in output")
    parser.add_argument("--S",   type=int, required=True,  help="Number of sequences to sample")
    parser.add_argument("--E",   type=float, required=True,
                        help="Edge density threshold [0.0, 1.0]")
    parser.add_argument("--seed", type=int, default=None,  help="Random seed for reproducibility")
    parser.add_argument("--out-fasta",         default="out_sequences.fasta",
                        help="Output FASTA path")
    parser.add_argument("--out-interactions",  default="out_interactions.txt",
                        help="Output interactions path")
    parser.add_argument("--out-motifs",        default="out_motifs.txt",
                        help="Output motifs path")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Basic validation
    if not (0.0 <= args.E <= 1.0):
        sys.exit("ERROR: E must be in [0.0, 1.0].")
    if args.d > args.l:
        sys.exit("ERROR: d cannot exceed l.")
    if args.K < 1:
        sys.exit("ERROR: K must be at least 1.")
    if args.S < 1:
        sys.exit("ERROR: S must be at least 1.")

    run(
        fasta_path=args.fasta,
        interactions_path=args.interactions,
        l=args.l,
        d=args.d,
        K=args.K,
        S=args.S,
        E=args.E,
        seed=args.seed,
        out_fasta=args.out_fasta,
        out_interactions=args.out_interactions,
        out_motifs=args.out_motifs,
    )
