#!/usr/bin/env python3
# =============================================================================
# align_fasta.py
#
# Align centroid FASTA records. The default backend is a small cross-platform
# Python center-star aligner; optional MAFFT/MUSCLE backends can be used when
# installed for more reliable multiple sequence alignment.
# =============================================================================

import argparse
import csv
import logging
import os
import shutil
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path
from time import perf_counter

from common import configure_logging
from config_schema import load_config_file
from fasta_io import count_fasta_records, parse_fasta, write_fasta

log = logging.getLogger(__name__)

_PAIRWISE_ALIGNER = None
ALIGNMENT_BACKENDS = {"python", "auto", "mafft", "muscle"}
WARN_ALIGNMENT_SEQUENCE_COUNT = 500
WARN_ALIGNMENT_BP = 2_000_000


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


def _first_available_backend():
    for backend in ["mafft", "muscle"]:
        if shutil.which(backend):
            return backend
    return None


def _run_and_log(cmd, stdout=None):
    log.info("Running: %s", " ".join(str(part) for part in cmd))
    if stdout is None:
        result = subprocess.run(cmd, capture_output=True, text=True)
    else:
        result = subprocess.run(cmd, stdout=stdout, stderr=subprocess.PIPE, text=True)
    if result.stdout:
        log.info(result.stdout.rstrip())
    if result.stderr:
        log.info(result.stderr.rstrip())
    return result


def _backend_executable(backend):
    if backend == "python":
        return sys.executable
    return shutil.which(backend) or ""


def _backend_version(backend):
    if backend == "python":
        return sys.version.split()[0]

    exe = shutil.which(backend)
    if not exe:
        return ""

    try:
        result = subprocess.run(
            [exe, "--version"], capture_output=True, text=True, timeout=10
        )
    except (OSError, subprocess.SubprocessError):
        return ""

    text = (result.stdout or result.stderr).strip().splitlines()
    return text[0] if text else ""


def write_alignment_metadata(path, row):
    if not path:
        return

    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "generated_at",
        "requested_backend",
        "backend_used",
        "fallback_used",
        "backend_executable",
        "backend_version",
        "n_input_sequences",
        "n_output_sequences",
        "input_total_bp",
        "elapsed_seconds",
    ]
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {"generated_at": datetime.now(UTC).isoformat(timespec="seconds"), **row}
        )


def run_mafft(input_path, output_path):
    exe = shutil.which("mafft")
    if not exe:
        raise FileNotFoundError("mafft executable not found")

    with open(output_path, "w", encoding="utf-8") as out_fh:
        result = _run_and_log([exe, "--auto", input_path], stdout=out_fh)
    return result.returncode == 0 and os.path.getsize(output_path) > 0


def run_muscle(input_path, output_path):
    exe = shutil.which("muscle")
    if not exe:
        raise FileNotFoundError("muscle executable not found")

    output_tmp = Path(output_path).with_suffix(Path(output_path).suffix + ".tmp")
    commands = [
        [exe, "-align", input_path, "-output", str(output_tmp)],
        [exe, "-in", input_path, "-out", str(output_tmp)],
    ]
    for cmd in commands:
        if output_tmp.exists():
            output_tmp.unlink()
        result = _run_and_log(cmd)
        if result.returncode == 0 and output_tmp.exists() and output_tmp.stat().st_size:
            os.replace(output_tmp, output_path)
            return True

    if output_tmp.exists():
        output_tmp.unlink()
    return False


def run_external_alignment(backend, input_path, output_path):
    runners = {"mafft": run_mafft, "muscle": run_muscle}
    return runners[backend](input_path, output_path)


def choose_backend(requested):
    if requested == "auto":
        backend = _first_available_backend()
        if backend:
            log.info("alignment_backend=auto selected %s", backend)
            return backend
        log.warning("alignment_backend=auto found no MAFFT/MUSCLE; using python")
        return "python"
    return requested


def align_records(records, input_path, output_path, requested_backend):
    backend = choose_backend(requested_backend)
    if backend == "python":
        aligned = center_star_align(records)
        write_fasta(aligned, output_path)
        return "python"

    try:
        ok = run_external_alignment(backend, input_path, output_path)
    except FileNotFoundError as exc:
        if requested_backend == "auto":
            log.warning("%s; falling back to python", exc)
            aligned = center_star_align(records)
            write_fasta(aligned, output_path)
            return "python"
        raise SystemExit(str(exc)) from exc

    if ok:
        return backend

    if requested_backend == "auto":
        log.warning("%s alignment failed; falling back to python", backend)
        aligned = center_star_align(records)
        write_fasta(aligned, output_path)
        return "python"

    raise SystemExit(f"{backend} alignment failed; see log for command output")


def parse_args():
    parser = argparse.ArgumentParser(description="Align FASTA records.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--config", help="Optional AmPrime config.yaml")
    parser.add_argument("--backend", choices=sorted(ALIGNMENT_BACKENDS))
    parser.add_argument("--metadata", help="Optional alignment metadata TSV")
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    configure_logging(args.log)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    cfg = load_config_file(args.config) if args.config else {}
    requested_backend = args.backend or cfg.get("alignment_backend", "python")
    if requested_backend not in ALIGNMENT_BACKENDS:
        raise SystemExit(
            "alignment backend must be one of: auto, mafft, muscle, python"
        )

    started = perf_counter()
    records = list(parse_fasta(args.input))
    n_in = len(records)
    if n_in < 2:
        shutil.copyfile(args.input, args.output)
        log.info("Only %d sequence(s); skipped alignment", n_in)
        elapsed = perf_counter() - started
        write_alignment_metadata(
            args.metadata,
            {
                "requested_backend": requested_backend,
                "backend_used": "skipped",
                "fallback_used": False,
                "backend_executable": "",
                "backend_version": "",
                "n_input_sequences": n_in,
                "n_output_sequences": n_in,
                "input_total_bp": sum(len(seq) for _, seq in records),
                "elapsed_seconds": round(elapsed, 3),
            },
        )
        return 0

    lengths = [len(seq) for _, seq in records]
    total_bp = sum(lengths)
    log.info(
        "Input size: %d sequence(s), %d bp total, length range %d-%d bp",
        n_in,
        total_bp,
        min(lengths),
        max(lengths),
    )
    if n_in > WARN_ALIGNMENT_SEQUENCE_COUNT or total_bp > WARN_ALIGNMENT_BP:
        log.warning(
            "Large alignment input. Python fallback alignment is intended for "
            "first-pass screening and may be slow or less accurate than a full MSA."
        )

    backend_used = align_records(records, args.input, args.output, requested_backend)

    n_out = count_fasta_records(args.output)
    elapsed = perf_counter() - started
    write_alignment_metadata(
        args.metadata,
        {
            "requested_backend": requested_backend,
            "backend_used": backend_used,
            "fallback_used": requested_backend == "auto" and backend_used == "python",
            "backend_executable": _backend_executable(backend_used),
            "backend_version": _backend_version(backend_used),
            "n_input_sequences": n_in,
            "n_output_sequences": n_out,
            "input_total_bp": total_bp,
            "elapsed_seconds": round(elapsed, 3),
        },
    )
    if backend_used == "python" and len({len(seq) for _, seq in records}) == 1:
        log.info(
            "All %d sequences have equal length; skipped pairwise alignment in %.2f s",
            n_out,
            elapsed,
        )
    else:
        log.info("Aligned %d sequences with %s in %.2f s", n_out, backend_used, elapsed)
    return 0


if __name__ == "__main__":
    sys.exit(main())
