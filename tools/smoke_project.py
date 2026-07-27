#!/usr/bin/env python3
"""Run fast project smoke checks that do not download data."""

import csv
import importlib.util
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "workflow" / "scripts"
SCRIPTS_PATH = str(SCRIPTS)

if SCRIPTS_PATH not in sys.path:
    sys.path.insert(0, SCRIPTS_PATH)


def smoke_env():
    env = os.environ.copy()
    pythonpath = env.get("PYTHONPATH")
    env["PYTHONPATH"] = (
        SCRIPTS_PATH if not pythonpath else os.pathsep.join([SCRIPTS_PATH, pythonpath])
    )
    return env


def load_script_module(module_name):
    module_path = SCRIPTS / f"{module_name}.py"
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {module_name} from {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def run_script_help(script_name):
    subprocess.run(
        [sys.executable, str(SCRIPTS / script_name), "--help"],
        cwd=ROOT,
        check=True,
        env=smoke_env(),
        stdout=subprocess.DEVNULL,
    )
    print(f"help ok {script_name}")


def check_config_validation():
    config_schema = load_script_module("config_schema")
    load_config_file = config_schema.load_config_file
    cfg = load_config_file(ROOT / "config" / "config.yaml")
    assert cfg["genus"]
    assert cfg["genes"]
    print("config validation ok")


def check_kmer_boundary():
    import numpy as np

    design_primers = load_script_module("design_primers")
    design_primers.np = np
    kmers = design_primers._build_kmers(
        np.array(list("acgt")),
        np.array(list("ACGT")),
        np.array([1, 1, 1, 1]),
        np.array([0.1, 0.2, 0.3, 0.4]),
        4,
    )
    assert len(kmers) == 1
    assert kmers[0]["degen"] == "ACGT"
    print("kmer boundary ok")


def check_primer_qc_cli():
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp = Path(tmp_dir)
        in_tsv = tmp / "primers.tsv"
        out_tsv = tmp / "checked.tsv"
        log = tmp / "qc.log"
        in_tsv.write_text(
            "primer_id\tfwd\trev\np1\tAAAAAA\tAAAAAA\np2\tAAAAAA\tTTTTTT\n",
            encoding="utf-8",
        )
        subprocess.run(
            [
                sys.executable,
                str(SCRIPTS / "check_primers.py"),
                "--in-tsv",
                str(in_tsv),
                "--out-tsv",
                str(out_tsv),
                "--max-heterodimer-dg",
                "-1",
                "--log",
                str(log),
            ],
            cwd=ROOT,
            check=True,
            env=smoke_env(),
        )
        with out_tsv.open(encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert [row["primer_id"] for row in rows] == ["p1"]
        assert rows[0]["heterodimer_dg"] == "0.0"
    print("primer QC cli ok")


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


def check_sequence_cli_steps():
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp = Path(tmp_dir)
        raw_fasta = tmp / "raw.fasta"
        centroids = tmp / "centroids.fasta"
        aligned = tmp / "aligned.fasta"
        cluster_log = tmp / "cluster.log"
        align_log = tmp / "align.log"

        raw_fasta.write_text(
            ">a\nACGTACGTACGT\n"
            ">b\nACGTACGTACGT\n"
            ">c\nACGTACGTTTGT\n",
            encoding="utf-8",
        )
        subprocess.run(
            [
                sys.executable,
                str(SCRIPTS / "cluster_fasta.py"),
                "--input",
                str(raw_fasta),
                "--output",
                str(centroids),
                "--identity",
                "0.97",
                "--log",
                str(cluster_log),
            ],
            cwd=ROOT,
            check=True,
            env=smoke_env(),
        )
        assert sum(1 for _ in parse_fasta(centroids)) == 2

        subprocess.run(
            [
                sys.executable,
                str(SCRIPTS / "align_fasta.py"),
                "--input",
                str(centroids),
                "--output",
                str(aligned),
                "--log",
                str(align_log),
            ],
            cwd=ROOT,
            check=True,
            env=smoke_env(),
        )
        aligned_lengths = {len(seq) for _, seq in parse_fasta(aligned)}
        assert len(aligned_lengths) == 1
    print("cluster and align cli ok")


def check_in_silico_pcr_cli():
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp = Path(tmp_dir)
        genome_dir = tmp / "genomes"
        genome_dir.mkdir()
        primers = tmp / "primers.tsv"
        out_tsv = tmp / "amplicons.tsv"
        log = tmp / "pcr.log"

        primers.write_text(
            "primer_id\tfwd\trev\np1\tATGC\tGCGT\n",
            encoding="utf-8",
        )
        (genome_dir / "genome.fna").write_text(
            ">contig1\nATGCGGGGGGGGGGACGC\n",
            encoding="utf-8",
        )
        subprocess.run(
            [
                sys.executable,
                str(SCRIPTS / "in_silico_pcr.py"),
                "--primers-tsv",
                str(primers),
                "--genome-dir",
                str(genome_dir),
                "--out-tsv",
                str(out_tsv),
                "--gene",
                "test",
                "--mismatch",
                "0",
                "--amplicon-min-len",
                "10",
                "--amplicon-max-len",
                "30",
                "--log",
                str(log),
            ],
            cwd=ROOT,
            check=True,
            env=smoke_env(),
        )
        with out_tsv.open(encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert rows[0]["n_genomes_amplified"] == "1"
    print("in silico PCR cli ok")


def main():
    for script_name in [
        "download_genomes.py",
        "extract_gene.py",
        "cluster_fasta.py",
        "align_fasta.py",
        "design_primers.py",
        "check_primers.py",
        "in_silico_pcr.py",
        "gene_report.py",
    ]:
        run_script_help(script_name)

    check_config_validation()
    check_kmer_boundary()
    check_primer_qc_cli()
    check_sequence_cli_steps()
    check_in_silico_pcr_cli()


if __name__ == "__main__":
    main()
