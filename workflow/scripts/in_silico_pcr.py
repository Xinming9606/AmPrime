#!/usr/bin/env python3
# =============================================================================
# in_silico_pcr.py
#
# Validate top-ranked primer pairs against full genome assemblies with a
# cross-platform Python primer scanner.
# =============================================================================

import argparse
import csv
import glob
import logging
import os
import sys
from bisect import bisect_left, bisect_right
from time import perf_counter

from config_schema import load_config_file
from fasta_io import parse_fasta

log = logging.getLogger(__name__)
WARN_GENOME_COUNT = 500
WARN_TOTAL_GENOME_BP = 500_000_000

OUT_COLS = [
    "validation_rank",
    "input_rank",
    "primer_id",
    "fwd",
    "rev",
    "n_genomes_amplified",
    "total_genomes",
    "amplification_rate",
    "mean_amplicon_len",
    "combined_score",
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
        writer = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


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
    sequence = sequence.upper()
    return amplicon_lengths_for_strands(
        (sequence, reverse_complement(sequence)), fwd, rev, mismatch, valid_lo, valid_hi
    )


def amplicon_lengths_for_strands(strands, fwd, rev, mismatch, valid_lo, valid_hi):
    rev_binding = reverse_complement(rev)
    len_fwd = len(fwd)
    len_rev = len(rev)
    lengths = []
    for strand in strands:
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
    parser.add_argument("--top-n", type=int)
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


def _float_or_default(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _make_candidate(row, input_rank):
    return {
        "input_rank": input_rank,
        "primer_id": row["primer_id"],
        "fwd": row["fwd"],
        "rev": row["rev"],
        "combined_score": _float_or_default(row.get("combined_score")),
        "amplified_genomes": 0,
        "amplicon_lengths": [],
    }


def _summarize_candidates(candidates, total_genomes):
    rows = []
    for candidate in candidates:
        lengths = candidate["amplicon_lengths"]
        rate = (
            candidate["amplified_genomes"] / total_genomes if total_genomes else 0.0
        )
        mean_len = sum(lengths) / len(lengths) if lengths else 0.0
        rows.append(
            {
                "validation_rank": 0,
                "input_rank": candidate["input_rank"],
                "primer_id": candidate["primer_id"],
                "fwd": candidate["fwd"],
                "rev": candidate["rev"],
                "n_genomes_amplified": candidate["amplified_genomes"],
                "total_genomes": total_genomes,
                "amplification_rate": round(rate, 4),
                "mean_amplicon_len": round(mean_len, 1),
                "combined_score": candidate["combined_score"],
            }
        )

    rows.sort(
        key=lambda row: (
            -_float_or_default(row["amplification_rate"]),
            -_float_or_default(row["combined_score"]),
            int(row["input_rank"]),
        )
    )
    for rank, row in enumerate(rows, 1):
        row["validation_rank"] = rank
    return rows


def main():
    args = parse_args()
    configure_logging(args.log)
    started = perf_counter()

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
    top_n = _required_param("pcr_top_n", _param(args.top_n, cfg, "pcr_top_n"))
    if top_n < 1:
        raise SystemExit("pcr_top_n must be a positive integer")

    len_margin = 100
    valid_lo = max(0, amplicon_min_len - len_margin)
    valid_hi = amplicon_max_len + len_margin

    log.info("Gene        : %s", args.gene)
    log.info("Primers TSV : %s", args.primers_tsv)
    log.info("Genome dir  : %s", args.genome_dir)
    log.info("Mismatch    : %d", mismatch)
    log.info("Top N       : %d", top_n)

    with open(args.primers_tsv, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        primer_rows = list(reader)

    if len(primer_rows) == 0:
        log.warning(
            "Primer TSV is empty - no primers to validate. Writing empty summary."
        )
        write_summary([], args.out_tsv)
        return 0

    candidates = [
        _make_candidate(row, input_rank)
        for input_rank, row in enumerate(primer_rows[:top_n], 1)
    ]
    log.info("Validating %d primer pair(s)", len(candidates))

    genomes = sorted(glob.glob(os.path.join(args.genome_dir, "*.fna")))
    total_genomes = len(genomes)
    log.info("Found %d genome files", total_genomes)
    if total_genomes > WARN_GENOME_COUNT:
        log.warning(
            "Large PCR input (%d genomes). Python scanning may be CPU-bound.",
            total_genomes,
        )

    if total_genomes == 0:
        log.error("No genome files found in %s", args.genome_dir)
        write_summary([], args.out_tsv)
        return 1

    total_bp_scanned = 0
    total_contigs = 0
    warned_total_bp = False
    for genome in genomes:
        genome_id = os.path.basename(genome)
        genome_bp = 0
        genome_contigs = 0
        per_candidate_hits = [[] for _ in candidates]
        for _, sequence in parse_fasta(genome):
            sequence = sequence.upper()
            genome_bp += len(sequence)
            genome_contigs += 1
            strands = (sequence, reverse_complement(sequence))
            for idx, candidate in enumerate(candidates):
                per_candidate_hits[idx].extend(
                    amplicon_lengths_for_strands(
                        strands,
                        candidate["fwd"],
                        candidate["rev"],
                        mismatch,
                        valid_lo,
                        valid_hi,
                    )
                )

        total_bp_scanned += genome_bp
        total_contigs += genome_contigs
        amplified_here = 0
        for candidate, lengths in zip(candidates, per_candidate_hits, strict=False):
            if not lengths:
                continue
            amplified_here += 1
            candidate["amplified_genomes"] += 1
            candidate["amplicon_lengths"].extend(lengths)

        if amplified_here:
            log.info(
                "  %s : amplified by %d pair(s); scanned %d contig(s), %d bp",
                genome_id,
                amplified_here,
                genome_contigs,
                genome_bp,
            )
        else:
            log.info(
                "  %s : no amplification; scanned %d contig(s), %d bp",
                genome_id,
                genome_contigs,
                genome_bp,
            )

        if not warned_total_bp and total_bp_scanned > WARN_TOTAL_GENOME_BP:
            log.warning(
                "PCR scan has exceeded %d bp. Consider stricter assembly_level "
                "or a smaller pcr_top_n for faster batch runs.",
                WARN_TOTAL_GENOME_BP,
            )
            warned_total_bp = True

    rows = _summarize_candidates(candidates, total_genomes)
    best = rows[0]
    log.info(
        "Best pair after PCR: %s amplified %s / %s genomes (rate %.3f)",
        best["primer_id"],
        best["n_genomes_amplified"],
        best["total_genomes"],
        _float_or_default(best["amplification_rate"]),
    )

    write_summary(rows, args.out_tsv)
    elapsed = perf_counter() - started
    log.info(
        "Scanned %d genome(s), %d contig(s), %d bp for %d candidate pair(s) in %.2f s",
        total_genomes,
        total_contigs,
        total_bp_scanned,
        len(candidates),
        elapsed,
    )
    log.info("Written to %s", args.out_tsv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
