#!/usr/bin/env python3
# =============================================================================
# align_fasta.py
#
# Align centroid FASTA records with MUSCLE. Inputs with fewer than two records
# are copied through unchanged so downstream reporting can continue.
# =============================================================================

import argparse
import logging
import os
import shutil
import subprocess
import sys

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
    n = 0
    with open(path, encoding="utf-8-sig") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def parse_args():
    parser = argparse.ArgumentParser(description="Align FASTA records with MUSCLE.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    configure_logging(args.log)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    n_in = count_fasta_records(args.input)
    if n_in < 2:
        shutil.copyfile(args.input, args.output)
        log.info("Only %d sequence(s); skipped MUSCLE alignment", n_in)
        return 0

    cmd = ["muscle", "-align", args.input, "-output", args.output]
    log.info("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.stdout:
        log.info(result.stdout.rstrip())
    if result.stderr:
        log.info(result.stderr.rstrip())
    if result.returncode != 0:
        log.error("MUSCLE failed with exit code %d", result.returncode)
        return result.returncode

    n_out = count_fasta_records(args.output)
    log.info("Aligned %d sequences", n_out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
