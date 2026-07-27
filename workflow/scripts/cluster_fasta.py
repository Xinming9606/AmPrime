#!/usr/bin/env python3
# =============================================================================
# cluster_fasta.py
#
# Dereplicate FASTA records with a small cross-platform Python implementation.
# Empty input produces an empty centroid FASTA so downstream per-gene reporting
# can continue.
# =============================================================================

import argparse
import logging
import os
import sys
from collections import defaultdict
from math import ceil, floor
from time import perf_counter

from fasta_io import count_fasta_records, parse_fasta, write_fasta

log = logging.getLogger(__name__)

_GLOBALXX_ALIGNER = None
WARN_SEQUENCE_COUNT = 1000
WARN_CENTROID_COUNT = 500


def configure_logging(log_path):
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def same_length_identity(seq_a, seq_b):
    matches = sum(a == b for a, b in zip(seq_a, seq_b, strict=False))
    return matches / len(seq_a)


def length_identity_upper_bound(len_a, len_b):
    if not len_a or not len_b:
        return 0.0
    return min(len_a, len_b) / max(len_a, len_b)


def globalxx_aligner():
    global _GLOBALXX_ALIGNER
    if _GLOBALXX_ALIGNER is None:
        from Bio.Align import PairwiseAligner

        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 1
        aligner.mismatch_score = 0
        aligner.open_gap_score = 0
        aligner.extend_gap_score = 0
        _GLOBALXX_ALIGNER = aligner
    return _GLOBALXX_ALIGNER


def sequence_identity(seq_a, seq_b):
    """Return an approximate global identity for clustering."""
    if not seq_a or not seq_b:
        return 0.0
    if len(seq_a) == len(seq_b):
        return same_length_identity(seq_a, seq_b)

    matches = globalxx_aligner().align(seq_a, seq_b).score
    return matches / max(len(seq_a), len(seq_b))


def cluster_records(records, identity):
    """Greedy centroid clustering in input order."""
    if not 0 < identity <= 1:
        raise ValueError("--identity must be greater than 0 and at most 1")

    centroids = []
    centroids_by_length = defaultdict(list)
    seen_sequences = set()

    for header, seq in records:
        seq_norm = seq.upper()
        if seq_norm in seen_sequences:
            continue

        seq_len = len(seq_norm)
        min_len = ceil(identity * seq_len)
        max_len = floor(seq_len / identity)
        redundant = False

        for length in range(min_len, max_len + 1):
            for centroid_seq in centroids_by_length.get(length, []):
                if (
                    length_identity_upper_bound(seq_len, len(centroid_seq)) >= identity
                    and sequence_identity(seq_norm, centroid_seq) >= identity
                ):
                    redundant = True
                    break
            if redundant:
                break

        if redundant:
            seen_sequences.add(seq_norm)
            continue

        centroids.append((header, seq))
        centroids_by_length[seq_len].append(seq_norm)
        seen_sequences.add(seq_norm)
    return centroids


def parse_args():
    parser = argparse.ArgumentParser(description="Cluster FASTA records.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--identity", required=True, type=float)
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    configure_logging(args.log)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    started = perf_counter()
    records = list(parse_fasta(args.input))
    n_in = len(records)
    if n_in == 0:
        open(args.output, "w", encoding="utf-8").close()
        log.info("No input sequences; wrote empty centroids FASTA")
        return 0

    total_bp = sum(len(seq) for _, seq in records)
    log.info("Input size: %d sequence(s), %d bp total", n_in, total_bp)
    if n_in > WARN_SEQUENCE_COUNT:
        log.warning(
            "Large cluster input (%d sequences). Greedy Python clustering may be slow.",
            n_in,
        )

    centroids = cluster_records(records, args.identity)
    if len(centroids) > WARN_CENTROID_COUNT:
        log.warning(
            "Many centroids retained (%d). Alignment and primer design may be slow.",
            len(centroids),
        )

    write_fasta(centroids, args.output)

    n_out = count_fasta_records(args.output)
    elapsed = perf_counter() - started
    log.info(
        "Clustered %d sequences into %d centroids in %.2f s", n_in, n_out, elapsed
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
