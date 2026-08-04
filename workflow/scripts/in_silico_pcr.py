#!/usr/bin/env python3
# =============================================================================
# in_silico_pcr.py
#
# Validate top-ranked primer pairs against full genome assemblies with SeqKit.
# =============================================================================

import argparse
import csv
import glob
import logging
import os
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter

from common import (
    config_param as _param,
    configure_logging,
    required_param as _required_param,
    reverse_complement,
)
from config_schema import load_config_file
from dependencies import ensure_tool
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

SPECIES_SUMMARY_COLS = ["metric", "value"]
SPECIES_TABLE_COLS = [
    "species",
    "n_amplicon_alleles",
    "shared_alleles",
    "has_inter_species_overlap",
]


def write_summary(rows, out_tsv):
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    with open(out_tsv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _canonical_sequence(sequence):
    reverse = reverse_complement(sequence)
    return min(sequence, reverse)


def _species_from_header(header, fallback):
    match = re.search(r"\[organism=([^\]]+)\]", header, re.IGNORECASE)
    if match:
        return " ".join(match.group(1).split())
    return fallback


def _write_seqkit_primer_file(path, candidates):
    """Write SeqKit's three-column primer-pair input file."""
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerows(
            (candidate["primer_id"], candidate["fwd"], candidate["rev"])
            for candidate in candidates
        )


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
    """Run SeqKit once for one genome and return all candidate observations."""
    executable, genome, primer_file, candidate_ids, mismatch, valid_lo, valid_hi = task
    genome_bp = 0
    genome_contigs = 0
    per_candidate_hits = [[] for _ in candidate_ids]
    per_candidate_species = [{} for _ in candidate_ids]
    genome_species = set()
    record_species = {}
    for header, sequence in parse_fasta(genome):
        species = _species_from_header(header, os.path.basename(genome))
        genome_species.add(species)
        record_id = header.lstrip(">\ufeff").split(maxsplit=1)[0]
        record_species[record_id] = species
        genome_bp += len(sequence)
        genome_contigs += 1
    command = [
        executable,
        "amplicon",
        "--primer-file",
        primer_file,
        "--max-mismatch",
        str(mismatch),
        "--bed",
        "--threads",
        "1",
        "--quiet",
        genome,
    ]
    completed = subprocess.run(  # noqa: S603 - executable came from PATH lookup.
        command,
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        details = (completed.stderr or completed.stdout).strip()
        raise RuntimeError(
            f"SeqKit amplicon failed for {genome} with exit code "
            f"{completed.returncode}: {details}"
        )

    candidate_index = {primer_id: idx for idx, primer_id in enumerate(candidate_ids)}
    for line in completed.stdout.splitlines():
        fields = line.split("\t", 6)
        if len(fields) < 7:
            continue
        record_id, start, end, primer_id, _, _, amplicon = fields
        idx = candidate_index.get(primer_id)
        if idx is None:
            continue
        try:
            length = int(end) - int(start)
        except ValueError:
            continue
        if not valid_lo <= length <= valid_hi:
            continue
        species = record_species.get(record_id, os.path.basename(genome))
        per_candidate_hits[idx].append(length)
        per_candidate_species[idx].setdefault(species, set()).add(
            _canonical_sequence(amplicon)
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
        for species, alleles in species_alleles.items():
            candidate["species_alleles"].setdefault(species, set()).update(alleles)
        if sum(len(alleles) for alleles in species_alleles.values()) > 1:
            candidate["multi_allele_genomes"] += 1

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
            "Large PCR input (%d genomes). SeqKit scanning may be CPU-bound.",
            total_genomes,
        )

    if total_genomes == 0:
        log.error("No genome files found in %s", args.genome_dir)
        write_summary([], args.out_tsv)
        _write_species_outputs(args.species_summary, args.species_tsv, None, 0)
        return 1

    executable = ensure_tool("seqkit")
    candidate_ids = [candidate["primer_id"] for candidate in candidates]
    workers = min(args.workers, total_genomes)
    log.info("Scanning genomes with %d worker process(es)", workers)
    with TemporaryDirectory(prefix="amprime-seqkit-") as temp_dir:
        primer_file = str(Path(temp_dir) / "primers.tsv")
        _write_seqkit_primer_file(primer_file, candidates)
        tasks = [
            (
                executable,
                genome,
                primer_file,
                candidate_ids,
                mismatch,
                valid_lo,
                valid_hi,
            )
            for genome in genomes
        ]
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
