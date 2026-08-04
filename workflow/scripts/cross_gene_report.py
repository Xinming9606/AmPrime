#!/usr/bin/env python3
"""Build a compact cross-gene comparison report from species summaries."""

from __future__ import annotations

import argparse
import csv
import os
from html import escape
from pathlib import Path


def _read_metrics(path: str) -> dict[str, str]:
    with open(path, newline="", encoding="utf-8") as fh:
        return {
            row["metric"]: row["value"]
            for row in csv.DictReader(fh, delimiter="\t")
            if row.get("metric")
        }


def _percent(metrics: dict[str, str], key: str) -> str:
    try:
        return f"{float(metrics.get(key, 0)) * 100:.1f}%"
    except ValueError:
        return "0.0%"


def _gene_from_path(path: str) -> str:
    name = Path(path).name
    suffix = "_species_summary.tsv"
    return name.removesuffix(suffix) if name.endswith(suffix) else Path(name).stem


def _build_html(summary_paths: list[str]) -> str:
    rows = []
    for path in summary_paths:
        metrics = _read_metrics(path)
        rows.append(
            {
                "gene": _gene_from_path(path),
                "amplified_genomes": metrics.get("amplified_genomes", "0"),
                "amplification_rate": _percent(metrics, "amplification_rate"),
                "amplified_species": metrics.get("amplified_species", "0"),
                "multi_allele_genomes": metrics.get("multi_allele_genomes", "0"),
                "overlap_species": metrics.get("overlap_species", "0"),
                "overlap_rate": _percent(metrics, "overlap_rate"),
                "unique_amplicon_alleles": metrics.get("unique_amplicon_alleles", "0"),
            }
        )
    rows.sort(key=lambda row: row["gene"])
    body_rows = "\n".join(
        "<tr>" + "".join(f"<td>{escape(row[key])}</td>" for key in row) + "</tr>"
        for row in rows
    )
    headers = [
        ("gene", "Gene"),
        ("amplified_genomes", "Amplified genomes"),
        ("amplification_rate", "Amplification rate"),
        ("amplified_species", "Amplified species"),
        ("multi_allele_genomes", "Multiple-allele genomes"),
        ("overlap_species", "Overlap species"),
        ("overlap_rate", "Overlap rate"),
        ("unique_amplicon_alleles", "Unique alleles"),
    ]
    header_html = "".join(f"<th>{escape(label)}</th>" for _, label in headers)
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>AmPrime cross-gene comparison</title>
<style>
body {{ font-family: system-ui, sans-serif; color: #24313a; margin: 2rem; }}
main {{ max-width: 1100px; margin: auto; }}
table {{ border-collapse: collapse; width: 100%; }}
th, td {{ border: 1px solid #d8dee3; padding: .55rem .7rem; text-align: left; }}
th {{ background: #eef3f5; }}
.note {{ color: #56636b; }}
</style>
</head>
<body><main>
<h1>AmPrime cross-gene comparison</h1>
<p class="note">Best validated pair for each configured gene.</p>
<table><thead><tr>{header_html}</tr></thead><tbody>{body_rows}</tbody></table>
</main></body>
</html>
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-tsv", nargs="+", required=True)
    parser.add_argument("--out-html", required=True)
    args = parser.parse_args()
    os.makedirs(os.path.dirname(args.out_html), exist_ok=True)
    Path(args.out_html).write_text(_build_html(args.summary_tsv), encoding="utf-8")


if __name__ == "__main__":
    main()
