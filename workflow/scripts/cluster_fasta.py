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
import warnings

log = logging.getLogger(__name__)


def configure_logging(log_path):
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def count_fasta_records(path):
    return sum(1 for _ in parse_fasta(path))


def parse_fasta(path):
    header = None
    seq_parts = []
    with open(path, encoding="utf-8-sig") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def write_fasta(records, path):
    with open(path, "w", encoding="utf-8") as fh:
        for header, seq in records:
            fh.write(header + "\n")
            fh.writelines(seq[i : i + 80] + "\n" for i in range(0, len(seq), 80))


def sequence_identity(seq_a, seq_b):
    """Return an approximate global identity for clustering."""
    if not seq_a or not seq_b:
        return 0.0

    if len(seq_a) == len(seq_b):
        matches = sum(
            a == b for a, b in zip(seq_a.upper(), seq_b.upper(), strict=False)
        )
        return matches / len(seq_a)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from Bio import pairwise2

    aln = pairwise2.align.globalxx(
        seq_a.upper(), seq_b.upper(), one_alignment_only=True
    )
    if not aln:
        return 0.0
    matches = aln[0].score
    return matches / max(len(seq_a), len(seq_b))


def cluster_records(records, identity):
    """Greedy centroid clustering in input order."""
    centroids = []
    for header, seq in records:
        if any(
            sequence_identity(seq, centroid_seq) >= identity
            for _, centroid_seq in centroids
        ):
            continue
        centroids.append((header, seq))
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

    n_in = count_fasta_records(args.input)
    if n_in == 0:
        open(args.output, "w", encoding="utf-8").close()
        log.info("No input sequences; wrote empty centroids FASTA")
        return 0

    records = list(parse_fasta(args.input))
    centroids = cluster_records(records, args.identity)
    write_fasta(centroids, args.output)

    n_out = count_fasta_records(args.output)
    log.info("Clustered %d sequences into %d centroids", n_in, n_out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
