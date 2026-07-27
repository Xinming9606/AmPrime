#!/usr/bin/env python3
# =============================================================================
# gene_report.py - self-contained HTML report for a single gene
#
# Snakemake script.  Requires markdown (conda-forge).
# Builds the body in Markdown, converts to HTML, and wraps in a styled page.
# =============================================================================

import base64
import csv
import logging
import os
from datetime import datetime
from html import escape

import markdown

# =============================================================================
# HTML page shell - loaded from separate .html file
# =============================================================================
_HTML_DIR = os.path.dirname(__file__)
with open(os.path.join(_HTML_DIR, "gene_report.html"), encoding="utf-8") as _fh:
    _PAGE = _fh.read()

# =============================================================================
# Markdown engine - tables turned on
# =============================================================================
_MD_EXTENSIONS = ["tables"]


# =============================================================================
# Helpers
# =============================================================================
def _read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return [dict(row) for row in csv.DictReader(fh, delimiter="\t")]


def _b64_png(path):
    with open(path, "rb") as fh:
        return base64.b64encode(fh.read()).decode()


def _num(row, key, default=0.0):
    value = row.get(key, default)
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _md_cell(value):
    return escape(str(value)).replace("|", r"\|")


def _md_code(value):
    return f"`{_md_cell(value)}`"


def _build_body(gene, genus, timestamp, primers, top_primer, pcr, diversity_img):
    """Build the report body as a Markdown string, using inline HTML only
    for styled elements (alerts, badges).  All control flow is plain Python."""

    out = []
    gene_text = _md_cell(gene)
    genus_text = _md_cell(genus)

    # ---- title -----------------------------------------------------------
    out.append(f"# {genus_text} - *{gene_text}* primer design")
    out.append("")
    out.append(f"*Generated {timestamp}*")
    out.append("")

    # ---- in silico PCR ---------------------------------------------------
    out.append("## In silico PCR validation")
    out.append("")
    if pcr:
        a = pcr
        rate = f"{_num(a, 'amplification_rate') * 100:.1f}%"
        out.append("|  |  |")
        out.append("|---|---|")
        out.append(f"| **Primer pair** | {_md_code(a.get('primer_id', ''))} |")
        out.append(f"| **Forward primer** | {_md_code(a.get('fwd', ''))} |")
        out.append(f"| **Reverse primer** | {_md_code(a.get('rev', ''))} |")
        out.append(
            f"| **Genomes amplified** | {_md_cell(a.get('n_genomes_amplified', ''))} / "
            f"{_md_cell(a.get('total_genomes', ''))}"
            f' <span class="badge badge-green">{rate}</span> |'
        )
        out.append(
            f"| **Mean amplicon length** | {_num(a, 'mean_amplicon_len'):.0f} bp |"
        )
    else:
        out.append('<div class="alert alert-info">')
        out.append(
            "<strong>No primer pair was validated.</strong> "
            "Either no primers passed the design filters, or the top pair "
            "amplified no genomes. See the candidate table below."
        )
        out.append("</div>")
    out.append("")

    # ---- recommended pair ------------------------------------------------
    out.append("## Recommended primer pair")
    out.append("")
    if top_primer:
        t = top_primer
        out.append("| Property | Forward | Reverse |")
        out.append("|---|---|---|")
        out.append(
            f"| Sequence | {_md_code(t.get('fwd', ''))} | {_md_code(t.get('rev', ''))} |"
        )
        out.append(
            f"| Position | {_md_cell(t.get('fwd_pos', ''))} | "
            f"{_md_cell(t.get('rev_pos', ''))} |"
        )
        out.append(
            f"| GC fraction | {_num(t, 'fwd_GC'):.2f} | {_num(t, 'rev_GC'):.2f} |"
        )
        out.append("")
        out.append(
            f"**Amplicon length:** {_md_cell(t.get('amplicon_len', ''))} bp &ensp;"
            f"**Pair diversity (Shannon):** {_num(t, 'pair_diversity'):.3f} &ensp;"
            f"**GC difference:** {_num(t, 'delta_GC'):.3f} &ensp;"
            f"**Combined score:** {_num(t, 'combined_score'):.3f}"
        )
    else:
        out.append('<div class="alert alert-warn">')
        out.append(
            "<strong>No candidate primers were found.</strong> "
            "Try raising <code>div_cut</code> or relaxing <code>GC_tol</code> "
            "in <code>config.yaml</code>."
        )
        out.append("</div>")
    out.append("")

    # ---- diversity plot --------------------------------------------------
    out.append("## Sequence diversity")
    out.append("")
    if diversity_img:
        out.append(
            "Per-position Shannon entropy across the alignment. "
            "Low-entropy (conserved) regions are good primer targets. "
            "Red rectangles mark the binding sites of the top-scoring primer pairs."
        )
        out.append("")
        out.append(
            f'<img src="data:image/png;base64,{diversity_img}" '
            f'alt="Diversity plot for {escape(str(gene), quote=True)}">'
        )
    else:
        out.append('<div class="alert alert-info">Diversity plot not available.</div>')
    out.append("")

    # ---- all candidates --------------------------------------------------
    out.append("## All candidate primer pairs")
    out.append("")
    if primers:
        keys = list(primers[0].keys())
        out.append("| " + " | ".join(_md_cell(key) for key in keys) + " |")
        out.append("|" + "|".join(["---"] * len(keys)) + "|")
        for row in primers:
            vals = [_md_cell(row.get(k, "")) for k in keys]
            out.append("| " + " | ".join(vals) + " |")
    else:
        out.append('<div class="alert alert-info">No candidates to display.</div>')
    out.append("")

    # ---- footer ----------------------------------------------------------
    out.append(f"*Primer pipeline - {genus_text}*")

    return "\n".join(out)


# =============================================================================
# Main
# =============================================================================
def main():
    log_path = snakemake.log[0]
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    log = logging.getLogger()

    gene = snakemake.params["gene"]
    genus = snakemake.params["genus"]

    # -- read inputs -------------------------------------------------------
    primers = (
        _read_tsv(snakemake.input["primers"])
        if os.path.isfile(snakemake.input["primers"])
        else []
    )
    amplicons = (
        _read_tsv(snakemake.input["amplicons"])
        if os.path.isfile(snakemake.input["amplicons"])
        else []
    )
    log.info("primers: %d rows, amplicons: %d rows", len(primers), len(amplicons))

    has_diversity = os.path.isfile(snakemake.input["diversity"])

    # -- build markdown body ------------------------------------------------
    body_md = _build_body(
        gene=gene,
        genus=genus,
        timestamp=datetime.now().strftime("%Y-%m-%d %H:%M"),
        primers=primers,
        top_primer=primers[0] if primers else None,
        pcr=amplicons[0] if amplicons else None,
        diversity_img=_b64_png(snakemake.input["diversity"]) if has_diversity else None,
    )

    # -- markdown to HTML, then wrap in page shell --------------------------
    body_html = markdown.markdown(body_md, extensions=_MD_EXTENSIONS)
    title = f"{genus} - {gene} primer design"
    html = _PAGE.replace("{TITLE}", escape(title)).replace("{CONTENT}", body_html)

    # -- write --------------------------------------------------------------
    out_html = snakemake.output[0]
    os.makedirs(os.path.dirname(out_html), exist_ok=True)
    with open(out_html, "w", encoding="utf-8") as fh:
        fh.write(html)

    log.info("Wrote report -> %s", out_html)


if __name__ == "__main__":
    main()
