#!/usr/bin/env python3
# =============================================================================
# extract_gene.py
#
# Extract sequences for a target gene from downloaded NCBI CDS and RNA FASTA
# files. Missing genes produce an empty FASTA so per-gene reports can still be
# generated in batch runs.
# =============================================================================

import argparse
import logging
import os
import re
from time import perf_counter

from common import configure_logging
from config_schema import load_config_file
from fasta_io import parse_fasta, write_fasta

log = logging.getLogger(__name__)
WARN_FASTA_FILE_COUNT = 1000


def header_matches(header, names):
    """Return True if a FASTA header matches any target gene name or alias."""
    h = header.lower()

    gene_tag = re.search(r"\[gene=([^\]]+)\]", h)
    if gene_tag and gene_tag.group(1).strip() in names:
        return True

    product_tag = re.search(r"\[product=([^\]]+)\]", h)
    if product_tag:
        product = product_tag.group(1).strip()
        if any(name in product for name in names):
            return True

    return False


def extract_from_dir(directory, names, label, gene):
    """Scan one FASTA directory and return matching (header, seq) records."""
    started = perf_counter()
    fna_files = sorted(
        os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".fna")
    )

    if not fna_files:
        log.warning("No .fna files found in %s dir: %s", label, directory)
        return []

    log.info("Scanning %d %s files ...", len(fna_files), label)
    if len(fna_files) > WARN_FASTA_FILE_COUNT:
        log.warning(
            "Large %s input (%d FASTA files). Gene extraction may be I/O-bound.",
            label,
            len(fna_files),
        )

    extracted = []
    n_records = 0
    n_bp = 0
    for fna in fna_files:
        genome_id = os.path.basename(fna)
        hits = []
        for hdr, seq in parse_fasta(fna):
            n_records += 1
            n_bp += len(seq)
            if header_matches(hdr, names):
                hits.append((hdr, seq))
        if hits:
            log.info("  %s : %d sequence(s) found", genome_id, len(hits))
            extracted.extend(hits)
        else:
            log.warning("  %s : gene '%s' not found - skipping", genome_id, gene)

    elapsed = perf_counter() - started
    log.info(
        "Scanned %d %s FASTA record(s), %d bp total, %d hit(s) in %.2f s",
        n_records,
        label,
        n_bp,
        len(extracted),
        elapsed,
    )
    return extracted


def parse_args():
    parser = argparse.ArgumentParser(description="Extract target gene sequences.")
    parser.add_argument("--cds-dir", required=True)
    parser.add_argument("--rna-dir", required=True)
    parser.add_argument("--out-fasta", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--config", help="Optional AmPrime config.yaml for aliases")
    parser.add_argument("--alias", action="append", default=[], dest="aliases")
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    configure_logging(args.log)

    aliases = list(args.aliases)
    if args.config:
        cfg = load_config_file(args.config)
        config_aliases = cfg.get("gene_aliases", {}).get(args.gene, [])
        aliases.extend(config_aliases)

    search_names = set([args.gene.lower()] + [a.lower() for a in aliases])

    log.info("Target gene : %s", args.gene)
    log.info("Aliases     : %s", aliases)
    log.info("Search names: %s", sorted(search_names))
    log.info("CDS dir     : %s", args.cds_dir)
    log.info("RNA dir     : %s", args.rna_dir)

    results = []
    if os.path.isdir(args.cds_dir):
        results.extend(extract_from_dir(args.cds_dir, search_names, "CDS", args.gene))
    else:
        log.warning("CDS directory does not exist: %s", args.cds_dir)

    if os.path.isdir(args.rna_dir):
        results.extend(extract_from_dir(args.rna_dir, search_names, "RNA", args.gene))
    else:
        log.warning("RNA directory does not exist: %s", args.rna_dir)

    log.info("Total sequences extracted: %d", len(results))
    if not results:
        log.warning(
            "No sequences found for gene '%s' in any genome. "
            "Writing an empty FASTA so downstream report generation can continue.",
            args.gene,
        )

    write_fasta(results, args.out_fasta)
    log.info("Written to %s", args.out_fasta)


if __name__ == "__main__":
    main()
