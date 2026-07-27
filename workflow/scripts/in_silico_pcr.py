#!/usr/bin/env python3
# =============================================================================
# in_silico_pcr.py
#
# Validate the top-ranked primer pair against full genome assemblies using
# `seqkit amplicon`.
# =============================================================================

import argparse
import csv
import glob
import logging
import os
import subprocess
import sys
import tempfile

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


def parse_args():
    parser = argparse.ArgumentParser(description="Validate primers with seqkit.")
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
        raise SystemExit(f"missing --{name.replace('_', '-')} or config setting: {name}")
    return value


def _param(cli_value, cfg, key):
    return cli_value if cli_value is not None else cfg.get(key)


def main():
    args = parse_args()
    configure_logging(args.log)

    cfg = load_config_file(args.config) if args.config else {}
    mismatch = _required_param("pcr_mismatch", _param(args.mismatch, cfg, "pcr_mismatch"))
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
        log.warning("Primer TSV is empty - no primers to validate. Writing empty summary.")
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

    primer_file = None
    try:
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as pf:
            pf.write(f"{primer_id}\t{fwd}\t{rev}\n")
            primer_file = pf.name
        log.info("Wrote seqkit primer file: %s", primer_file)

        amplified_genomes = 0
        amplicon_lengths = []

        for genome in genomes:
            genome_id = os.path.basename(genome)
            cmd = [
                "seqkit",
                "amplicon",
                "-p",
                primer_file,
                "-m",
                str(mismatch),
                "--bed",
                genome,
            ]
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            except subprocess.CalledProcessError as exc:
                log.warning(
                    "  %s : seqkit failed (%s) - skipping",
                    genome_id,
                    exc.stderr.strip(),
                )
                continue

            lines = [ln for ln in result.stdout.splitlines() if ln.strip()]
            valid_here = []
            for line in lines:
                fields = line.split("\t")
                if len(fields) < 3:
                    continue
                try:
                    amp_len = int(fields[2]) - int(fields[1])
                except ValueError:
                    continue
                if valid_lo <= amp_len <= valid_hi:
                    valid_here.append(amp_len)

            if valid_here:
                amplified_genomes += 1
                amplicon_lengths.extend(valid_here)
                log.info(
                    "  %s : amplified (%d valid product[s], %d raw)",
                    genome_id,
                    len(valid_here),
                    len(lines),
                )
            elif lines:
                log.info(
                    "  %s : only spurious products (%d raw, none in %d-%d bp) - not counted",
                    genome_id,
                    len(lines),
                    valid_lo,
                    valid_hi,
                )
            else:
                log.info("  %s : no amplification", genome_id)
    finally:
        if primer_file and os.path.exists(primer_file):
            os.unlink(primer_file)

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
