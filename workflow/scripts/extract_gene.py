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

log = logging.getLogger(__name__)


def configure_logging(log_path):
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def parse_fasta(filepath):
    """Yield (header, sequence) tuples from a FASTA file."""
    header = None
    seq_parts = []
    with open(filepath, encoding="utf-8") as fh:
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
    fna_files = sorted(
        os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".fna")
    )

    if not fna_files:
        log.warning("No .fna files found in %s dir: %s", label, directory)
        return []

    log.info("Scanning %d %s files ...", len(fna_files), label)

    extracted = []
    for fna in fna_files:
        genome_id = os.path.basename(fna)
        hits = [
            (hdr, seq) for hdr, seq in parse_fasta(fna) if header_matches(hdr, names)
        ]
        if hits:
            log.info("  %s : %d sequence(s) found", genome_id, len(hits))
            extracted.extend(hits)
        else:
            log.warning("  %s : gene '%s' not found - skipping", genome_id, gene)

    return extracted


def write_fasta(records, out_fasta):
    os.makedirs(os.path.dirname(out_fasta), exist_ok=True)
    with open(out_fasta, "w", encoding="utf-8") as fh:
        for header, seq in records:
            fh.write(header + "\n")
            fh.writelines(seq[i : i + 80] + "\n" for i in range(0, len(seq), 80))


def parse_args():
    parser = argparse.ArgumentParser(description="Extract target gene sequences.")
    parser.add_argument("--cds-dir", required=True)
    parser.add_argument("--rna-dir", required=True)
    parser.add_argument("--out-fasta", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--alias", action="append", default=[], dest="aliases")
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    configure_logging(args.log)

    search_names = set([args.gene.lower()] + [a.lower() for a in args.aliases])

    log.info("Target gene : %s", args.gene)
    log.info("Aliases     : %s", args.aliases)
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
