#!/usr/bin/env python3
"""Download the functional-test genomes and package them as a tar archive."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOWNLOAD_SCRIPT = ROOT / "workflow" / "scripts" / "download_genomes.py"
DEFAULT_CONFIG = ROOT / "config" / "config.yaml"
DEFAULT_OUTPUT = ROOT / "data" / "borrelia-genomes.tar.gz"


def download_archive(output: Path, genus: str, config: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    partial_output = output.with_name(output.name + ".partial")
    if partial_output.exists():
        partial_output.unlink()

    with tempfile.TemporaryDirectory(prefix="amprime-test-data-") as temp_dir:
        root = Path(temp_dir)
        genomes = root / "genomes"
        log_path = root / "download.log"
        command = [
            sys.executable,
            str(DOWNLOAD_SCRIPT),
            "--config",
            str(config),
            "--genus",
            genus,
            "--genomic-dir",
            str(genomes / "genomic"),
            "--cds-dir",
            str(genomes / "cds"),
            "--rna-dir",
            str(genomes / "rna"),
            "--manifest",
            str(genomes / "download_manifest.tsv"),
            "--log",
            str(log_path),
        ]
        subprocess.run(command, cwd=ROOT, check=True, env=_project_env())

        with tarfile.open(partial_output, "w:gz") as archive:
            archive.add(genomes, arcname="genomes")
    os.replace(partial_output, output)


def _project_env() -> dict[str, str]:
    env = os.environ.copy()
    scripts_dir = str(DOWNLOAD_SCRIPT.parent)
    pythonpath = env.get("PYTHONPATH")
    env["PYTHONPATH"] = (
        scripts_dir if not pythonpath else f"{scripts_dir}{os.pathsep}{pythonpath}"
    )
    return env


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genus", default="Borrelia")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    download_archive(args.output, args.genus, args.config)
    print(f"downloaded test dataset: {args.output}")


if __name__ == "__main__":
    main()
