#!/usr/bin/env python3
# =============================================================================
# align_fasta.py
#
# Align centroid FASTA records with a small cross-platform Python center-star
# aligner. Inputs with fewer than two records are copied through unchanged so
# downstream reporting can continue.
# =============================================================================

import argparse
import logging
import os
import shutil
import sys

log = logging.getLogger(__name__)

_PAIRWISE_ALIGNER = None


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


def pairwise_aligner():
    global _PAIRWISE_ALIGNER
    if _PAIRWISE_ALIGNER is None:
        from Bio.Align import PairwiseAligner

        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -0.5
        _PAIRWISE_ALIGNER = aligner
    return _PAIRWISE_ALIGNER


def pairwise_align(reference, sequence):
    alignments = pairwise_aligner().align(reference.upper(), sequence.upper())
    if not alignments:
        return reference.upper(), sequence.upper()
    alignment = alignments[0]
    return str(alignment[0]), str(alignment[1])


def merge_alignment_rows(ref_row, rows, new_ref_row, new_seq_row):
    merged_ref = []
    merged_rows = [[] for _ in rows]
    merged_new = []
    i = 0
    j = 0

    while i < len(ref_row) or j < len(new_ref_row):
        if i < len(ref_row) and ref_row[i] == "-":
            merged_ref.append("-")
            for idx, row in enumerate(rows):
                merged_rows[idx].append(row[i])
            merged_new.append("-")
            i += 1
        elif j < len(new_ref_row) and new_ref_row[j] == "-":
            merged_ref.append("-")
            for row in merged_rows:
                row.append("-")
            merged_new.append(new_seq_row[j])
            j += 1
        else:
            merged_ref.append(ref_row[i])
            for idx, row in enumerate(rows):
                merged_rows[idx].append(row[i])
            merged_new.append(new_seq_row[j])
            i += 1
            j += 1

    return (
        "".join(merged_ref),
        ["".join(row) for row in merged_rows],
        "".join(merged_new),
    )


def center_star_align(records):
    if len(records) < 2:
        return records

    if len({len(seq) for _, seq in records}) == 1:
        return [(header, seq.upper()) for header, seq in records]

    ref_index = max(range(len(records)), key=lambda idx: len(records[idx][1]))
    ref_header, ref_seq = records[ref_index]
    ordered = [records[ref_index]] + [
        record for idx, record in enumerate(records) if idx != ref_index
    ]

    ref_row = ref_seq.upper()
    rows = [ref_row]
    headers = [ref_header]

    for header, seq in ordered[1:]:
        new_ref_row, new_seq_row = pairwise_align(ref_seq, seq)
        ref_row, rows, new_seq_row = merge_alignment_rows(
            ref_row, rows, new_ref_row, new_seq_row
        )
        rows[0] = ref_row
        rows.append(new_seq_row)
        headers.append(header)

    return list(zip(headers, rows, strict=False))


def parse_args():
    parser = argparse.ArgumentParser(description="Align FASTA records.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    configure_logging(args.log)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    records = list(parse_fasta(args.input))
    n_in = len(records)
    if n_in < 2:
        shutil.copyfile(args.input, args.output)
        log.info("Only %d sequence(s); skipped alignment", n_in)
        return 0

    aligned = center_star_align(records)
    write_fasta(aligned, args.output)

    n_out = count_fasta_records(args.output)
    if len({len(seq) for _, seq in records}) == 1:
        log.info(
            "All %d sequences have equal length; skipped pairwise alignment", n_out
        )
    else:
        log.info("Aligned %d sequences", n_out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
