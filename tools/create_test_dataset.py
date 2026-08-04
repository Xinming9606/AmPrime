#!/usr/bin/env python3
"""Create a tiny deterministic FASTA archive for offline CI functional tests."""

from __future__ import annotations

import argparse
import tarfile
import tempfile
from pathlib import Path


def create_archive(output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as temp_dir:
        genomes = Path(temp_dir) / "genomes"
        contents = {
            "genomic": ">synthetic_genome\n" + "ACGT" * 100 + "\n",
            "cds": (
                ">synthetic_cds_1 [gene=recG]\n"
                + "ACGT" * 10
                + "\n>synthetic_cds_2 [gene=recG]\n"
                + "ACGA" * 10
                + "\n"
            ),
            "rna": ">synthetic_rna [gene=other]\n" + "ACGT" * 30 + "\n",
        }
        for label, content in contents.items():
            directory = genomes / label
            directory.mkdir(parents=True)
            (directory / f"synthetic_{label}.fna").write_text(content, encoding="utf-8")

        with tarfile.open(output, "w:gz") as archive:
            archive.add(genomes, arcname="genomes")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data") / "ci-borrelia-genomes.tar.gz",
    )
    args = parser.parse_args()
    create_archive(args.output)
    print(f"created test dataset: {args.output}")


if __name__ == "__main__":
    main()
