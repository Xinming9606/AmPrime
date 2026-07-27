#!/usr/bin/env python3
# =============================================================================
# in_silico_pcr.py
#
# Validate the top-ranked primer pair against full genome assemblies with a
# cross-platform Python primer scanner.
# =============================================================================

import argparse
import csv
import glob
import logging
import os
import sys
from bisect import bisect_left, bisect_right

from config_schema import load_config_file

log = logging.getLogger(__name__)

OUT_COLS = [
    "primer_id",
    "fwd",
    "rev",
    "n_genomes_amplified",
    "total_genomes",
    "amplification_rate",
    "mean_amplicon_len",
]

IUPAC_BASES = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "M": frozenset("AC"),
    "K": frozenset("GT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}

IUPAC_COMPLEMENT = str.maketrans(
    "ACGTRYMKSWHBVDNacgtrymkswhbvdn", "TGCAYRKMSWDVBHNtgcayrkmswdvbhn"
)


def configure_logging(log_path):
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def write_summary(rows, out_tsv):
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    with open(out_tsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(OUT_COLS)
        for row in rows:
            w.writerow(row)


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


def reverse_complement(seq):
    return seq.translate(IUPAC_COMPLEMENT)[::-1]


def base_matches(primer_base, target_base):
    primer_set = IUPAC_BASES.get(primer_base.upper(), frozenset())
    target_set = IUPAC_BASES.get(target_base.upper(), frozenset())
    return bool(primer_set & target_set)


def primer_window_mismatches(primer, sequence, start, max_mismatch):
    mismatches = 0
    for offset, primer_base in enumerate(primer):
        target_base = sequence[start + offset]
        if base_matches(primer_base, target_base):
            continue
        mismatches += 1
        if mismatches > max_mismatch:
            return mismatches
    return mismatches


def find_primer_sites(sequence, primer, max_mismatch):
    primer = primer.upper()
    sequence = sequence.upper()
    k = len(primer)
    for start in range(len(sequence) - k + 1):
        if (
            primer_window_mismatches(primer, sequence, start, max_mismatch)
            <= max_mismatch
        ):
            yield start


def amplicon_lengths_for_sequence(sequence, fwd, rev, mismatch, valid_lo, valid_hi):
    rev_binding = reverse_complement(rev)
    len_fwd = len(fwd)
    len_rev = len(rev)
    lengths = []
    for strand in [sequence, reverse_complement(sequence)]:
        fwd_sites = list(find_primer_sites(strand, fwd, mismatch))
        rev_sites = list(find_primer_sites(strand, rev_binding, mismatch))
        for fwd_start in fwd_sites:
            lo = max(fwd_start + len_fwd, fwd_start + valid_lo - len_rev)
            hi = fwd_start + valid_hi - len_rev
            start = bisect_left(rev_sites, lo)
            end = bisect_right(rev_sites, hi)
            lengths.extend(
                rev_start + len_rev - fwd_start for rev_start in rev_sites[start:end]
            )
    return lengths


def parse_args():
    parser = argparse.ArgumentParser(description="Validate primers in silico.")
    parser.add_argument("--primers-tsv", required=True)
    parser.add_argument("--genome-dir", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--config", help="Optional AmPrime config.yaml")
    parser.add_argument("--mismatch", type=int)
    parser.add_argument("--amplicon-min-len", type=int)
    parser.add_argument("--amplicon-max-len", type=int)
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def _required_param(name, value):
    if value is None:
        raise SystemExit(
            f"missing --{name.replace('_', '-')} or config setting: {name}"
        )
    return value


def _param(cli_value, cfg, key):
    return cli_value if cli_value is not None else cfg.get(key)


def main():
    args = parse_args()
    configure_logging(args.log)

    cfg = load_config_file(args.config) if args.config else {}
    mismatch = _required_param(
        "pcr_mismatch", _param(args.mismatch, cfg, "pcr_mismatch")
    )
    amplicon_min_len = _required_param(
        "amplicon_min_len", _param(args.amplicon_min_len, cfg, "amplicon_min_len")
    )
    amplicon_max_len = _required_param(
        "amplicon_max_len", _param(args.amplicon_max_len, cfg, "amplicon_max_len")
    )

    len_margin = 100
    valid_lo = max(0, amplicon_min_len - len_margin)
    valid_hi = amplicon_max_len + len_margin

    log.info("Gene        : %s", args.gene)
    log.info("Primers TSV : %s", args.primers_tsv)
    log.info("Genome dir  : %s", args.genome_dir)
    log.info("Mismatch    : %d", mismatch)

    with open(args.primers_tsv, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        primer_rows = list(reader)

    if len(primer_rows) == 0:
        log.warning(
            "Primer TSV is empty - no primers to validate. Writing empty summary."
        )
        write_summary([], args.out_tsv)
        return 0

    top = primer_rows[0]
    primer_id = top["primer_id"]
    fwd = top["fwd"]
    rev = top["rev"]
    log.info("Top primer pair: %s  fwd=%s  rev=%s", primer_id, fwd, rev)

    genomes = sorted(glob.glob(os.path.join(args.genome_dir, "*.fna")))
    total_genomes = len(genomes)
    log.info("Found %d genome files", total_genomes)

    if total_genomes == 0:
        log.error("No genome files found in %s", args.genome_dir)
        write_summary([], args.out_tsv)
        return 1

    amplified_genomes = 0
    amplicon_lengths = []

    for genome in genomes:
        genome_id = os.path.basename(genome)
        valid_here = []
        for _, sequence in parse_fasta(genome):
            valid_here.extend(
                amplicon_lengths_for_sequence(
                    sequence, fwd, rev, mismatch, valid_lo, valid_hi
                )
            )

        if valid_here:
            amplified_genomes += 1
            amplicon_lengths.extend(valid_here)
            log.info(
                "  %s : amplified (%d valid product[s])", genome_id, len(valid_here)
            )
        else:
            log.info("  %s : no amplification", genome_id)

    rate = amplified_genomes / total_genomes if total_genomes else 0.0
    mean_len = sum(amplicon_lengths) / len(amplicon_lengths) if amplicon_lengths else 0

    log.info(
        "Amplified %d / %d genomes (rate %.3f), mean amplicon length %.1f bp",
        amplified_genomes,
        total_genomes,
        rate,
        mean_len,
    )

    write_summary(
        [
            [
                primer_id,
                fwd,
                rev,
                amplified_genomes,
                total_genomes,
                round(rate, 4),
                round(mean_len, 1),
            ]
        ],
        args.out_tsv,
    )
    log.info("Written to %s", args.out_tsv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
