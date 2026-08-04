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
import re
import sys
from bisect import bisect_left, bisect_right
from concurrent.futures import ProcessPoolExecutor
from time import perf_counter

from common import (
    config_param as _param,
    configure_logging,
    required_param as _required_param,
    reverse_complement,
)
from config_schema import load_config_file
from fasta_io import parse_fasta

log = logging.getLogger(__name__)
WARN_GENOME_COUNT = 500
WARN_TOTAL_GENOME_BP = 500_000_000
MIN_SEED_LEN = 4
UNAMBIGUOUS_BASES = frozenset("ACGT")

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

SPECIES_SUMMARY_COLS = ["metric", "value"]
SPECIES_TABLE_COLS = [
    "species",
    "n_amplicon_alleles",
    "shared_alleles",
    "has_inter_species_overlap",
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


def write_summary(rows, out_tsv):
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    with open(out_tsv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


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


def _only_unambiguous_dna(seq):
    return all(base in UNAMBIGUOUS_BASES for base in seq)


def _seed_segments(primer, max_mismatch):
    """Split primer into max_mismatch + 1 exact-match seeds.

    Pigeonhole filtering is only safe when every segment can be searched
    exactly. If any segment is too short or contains an ambiguous base, callers
    should fall back to the full sliding-window scan.
    """

    if max_mismatch < 0 or max_mismatch >= len(primer):
        return None

    n_segments = max_mismatch + 1
    base_len, remainder = divmod(len(primer), n_segments)
    segments = []
    offset = 0
    for idx in range(n_segments):
        seg_len = base_len + (1 if idx < remainder else 0)
        seed = primer[offset : offset + seg_len]
        if len(seed) < MIN_SEED_LEN or not _only_unambiguous_dna(seed):
            return None
        segments.append((offset, seed))
        offset += seg_len
    return segments


def _seed_candidate_starts(sequence, primer, max_mismatch):
    if len(sequence) < len(primer):
        return []
    if not _only_unambiguous_dna(sequence):
        return None

    segments = _seed_segments(primer, max_mismatch)
    if segments is None:
        return None

    max_start = len(sequence) - len(primer)
    starts = set()
    for offset, seed in segments:
        found_at = sequence.find(seed)
        while found_at != -1:
            window_start = found_at - offset
            if 0 <= window_start <= max_start:
                starts.add(window_start)
            found_at = sequence.find(seed, found_at + 1)
    return sorted(starts)


def find_primer_sites(sequence, primer, max_mismatch):
    primer = primer.upper()
    sequence = sequence.upper()
    k = len(primer)
    candidate_starts = _seed_candidate_starts(sequence, primer, max_mismatch)
    starts = (
        candidate_starts
        if candidate_starts is not None
        else range(len(sequence) - k + 1)
    )
    for start in starts:
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


def amplicon_sequences_for_strands(strands, fwd, rev, mismatch, valid_lo, valid_hi):
    """Return primer-oriented amplicon sequences for both DNA strands."""
    rev_binding = reverse_complement(rev)
    len_fwd = len(fwd)
    len_rev = len(rev)
    sequences = []
    for strand in strands:
        fwd_sites = list(find_primer_sites(strand, fwd, mismatch))
        rev_sites = list(find_primer_sites(strand, rev_binding, mismatch))
        for fwd_start in fwd_sites:
            lo = max(fwd_start + len_fwd, fwd_start + valid_lo - len_rev)
            hi = fwd_start + valid_hi - len_rev
            start = bisect_left(rev_sites, lo)
            end = bisect_right(rev_sites, hi)
            sequences.extend(
                strand[fwd_start : rev_start + len_rev]
                for rev_start in rev_sites[start:end]
            )
    return sequences


def _canonical_sequence(sequence):
    reverse = reverse_complement(sequence)
    return min(sequence, reverse)


def _species_from_header(header, fallback):
    match = re.search(r"\[organism=([^\]]+)\]", header, re.IGNORECASE)
    if match:
        return " ".join(match.group(1).split())
    return fallback


def _write_species_outputs(summary_path, species_path, candidate, total_genomes):
    summary = candidate.get("species_alleles", {}) if candidate else {}
    all_species = candidate.get("all_species", set()) if candidate else set()
    amplified_genomes = candidate.get("amplified_genomes", 0) if candidate else 0
    multi_allele_genomes = candidate.get("multi_allele_genomes", 0) if candidate else 0
    shared_alleles = {
        allele
        for alleles in summary.values()
        for allele in alleles
        if sum(allele in values for values in summary.values()) > 1
    }
    overlap_species = sum(
        any(allele in shared_alleles for allele in alleles)
        for alleles in summary.values()
    )
    amplified_species = len(summary)
    metrics = {
        "total_genomes": total_genomes,
        "amplified_genomes": amplified_genomes,
        "amplification_rate": (
            round(amplified_genomes / total_genomes, 4) if total_genomes else 0.0
        ),
        "total_species": len(all_species),
        "amplified_species": amplified_species,
        "multi_allele_genomes": multi_allele_genomes,
        "multi_allele_rate": (
            round(multi_allele_genomes / amplified_genomes, 4)
            if amplified_genomes
            else 0.0
        ),
        "overlap_species": overlap_species,
        "overlap_rate": (
            round(overlap_species / amplified_species, 4) if amplified_species else 0.0
        ),
        "unique_amplicon_alleles": len(
            {allele for alleles in summary.values() for allele in alleles}
        ),
    }
    if summary_path:
        os.makedirs(os.path.dirname(summary_path), exist_ok=True)
        with open(summary_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=SPECIES_SUMMARY_COLS, delimiter="\t")
            writer.writeheader()
            writer.writerows(
                {"metric": key, "value": value} for key, value in metrics.items()
            )
    if species_path:
        os.makedirs(os.path.dirname(species_path), exist_ok=True)
        with open(species_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=SPECIES_TABLE_COLS, delimiter="\t")
            writer.writeheader()
            for species, alleles in sorted(summary.items()):
                shared = sum(allele in shared_alleles for allele in alleles)
                writer.writerow(
                    {
                        "species": species,
                        "n_amplicon_alleles": len(alleles),
                        "shared_alleles": shared,
                        "has_inter_species_overlap": bool(shared),
                    }
                )


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
    parser.add_argument(
        "--workers", type=int, default=1, help="Genome-scanning worker processes."
    )
    parser.add_argument("--species-summary")
    parser.add_argument("--species-tsv")
    parser.add_argument("--log", required=True)
    return parser.parse_args()


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
        "species_alleles": {},
        "all_species": set(),
        "multi_allele_genomes": 0,
    }


def _summarize_candidates(candidates, total_genomes):
    rows = []
    for candidate in candidates:
        lengths = candidate["amplicon_lengths"]
        rate = candidate["amplified_genomes"] / total_genomes if total_genomes else 0.0
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


def _scan_genome(task):
    """Scan one genome in a worker process and return raw candidate hits."""
    genome, candidate_specs, mismatch, valid_lo, valid_hi = task
    genome_bp = 0
    genome_contigs = 0
    per_candidate_hits = [[] for _ in candidate_specs]
    per_candidate_species = [{} for _ in candidate_specs]
    genome_species = set()
    for header, sequence in parse_fasta(genome):
        species = _species_from_header(header, os.path.basename(genome))
        genome_species.add(species)
        sequence = sequence.upper()
        genome_bp += len(sequence)
        genome_contigs += 1
        strands = (sequence, reverse_complement(sequence))
        for idx, (fwd, rev) in enumerate(candidate_specs):
            amplicons = amplicon_sequences_for_strands(
                strands, fwd, rev, mismatch, valid_lo, valid_hi
            )
            if not amplicons:
                continue
            per_candidate_hits[idx].extend(len(amplicon) for amplicon in amplicons)
            species_alleles = per_candidate_species[idx].setdefault(species, set())
            species_alleles.update(
                _canonical_sequence(amplicon) for amplicon in amplicons
            )
    return (
        os.path.basename(genome),
        genome_bp,
        genome_contigs,
        per_candidate_hits,
        per_candidate_species,
        genome_species,
    )


def _record_genome_result(scan_result, candidates):
    (
        genome_id,
        genome_bp,
        genome_contigs,
        per_candidate_hits,
        per_candidate_species,
        genome_species,
    ) = scan_result
    amplified_here = 0
    for candidate, lengths, species_alleles in zip(
        candidates, per_candidate_hits, per_candidate_species, strict=False
    ):
        candidate["all_species"].update(genome_species)
        if not lengths:
            continue
        amplified_here += 1
        candidate["amplified_genomes"] += 1
        candidate["amplicon_lengths"].extend(lengths)
        if sum(len(alleles) for alleles in species_alleles.values()) > 1:
            candidate["multi_allele_genomes"] += 1
        for species, alleles in species_alleles.items():
            candidate["species_alleles"].setdefault(species, set()).update(alleles)

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
    return genome_bp, genome_contigs


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
    if args.workers < 1:
        raise SystemExit("workers must be a positive integer")

    len_margin = 100
    valid_lo = max(0, amplicon_min_len - len_margin)
    valid_hi = amplicon_max_len + len_margin

    log.info("Gene        : %s", args.gene)
    log.info("Primers TSV : %s", args.primers_tsv)
    log.info("Genome dir  : %s", args.genome_dir)
    log.info("Mismatch    : %d", mismatch)
    log.info("Top N       : %d", top_n)
    log.info("Workers     : %d", args.workers)

    with open(args.primers_tsv, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        primer_rows = list(reader)

    if len(primer_rows) == 0:
        log.warning(
            "Primer TSV is empty - no primers to validate. Writing empty summary."
        )
        write_summary([], args.out_tsv)
        _write_species_outputs(args.species_summary, args.species_tsv, None, 0)
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
        _write_species_outputs(args.species_summary, args.species_tsv, None, 0)
        return 1

    candidate_specs = [(candidate["fwd"], candidate["rev"]) for candidate in candidates]
    tasks = [
        (genome, candidate_specs, mismatch, valid_lo, valid_hi) for genome in genomes
    ]
    workers = min(args.workers, total_genomes)
    log.info("Scanning genomes with %d worker process(es)", workers)
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            scan_results = executor.map(_scan_genome, tasks)
            totals = [
                _record_genome_result(result, candidates) for result in scan_results
            ]
    else:
        totals = [
            _record_genome_result(_scan_genome(task), candidates) for task in tasks
        ]

    total_bp_scanned = sum(bp for bp, _ in totals)
    total_contigs = sum(contigs for _, contigs in totals)
    if total_bp_scanned > WARN_TOTAL_GENOME_BP:
        log.warning(
            "PCR scan has exceeded %d bp. Consider stricter assembly_level "
            "or a smaller pcr_top_n for faster batch runs.",
            WARN_TOTAL_GENOME_BP,
        )

    rows = _summarize_candidates(candidates, total_genomes)
    best = rows[0]
    best_candidate = next(
        candidate
        for candidate in candidates
        if candidate["primer_id"] == best["primer_id"]
    )
    log.info(
        "Best pair after PCR: %s amplified %s / %s genomes (rate %.3f)",
        best["primer_id"],
        best["n_genomes_amplified"],
        best["total_genomes"],
        _float_or_default(best["amplification_rate"]),
    )

    write_summary(rows, args.out_tsv)
    _write_species_outputs(
        args.species_summary,
        args.species_tsv,
        best_candidate,
        total_genomes,
    )
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
