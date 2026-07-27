#!/usr/bin/env python3
# =============================================================================
# design_primers.py - pure Python primer design implementation
#
# Reads a multiple sequence alignment (FASTA), computes per-position Shannon
# entropy, builds degenerate (IUPAC) consensus kmer windows, evaluates all
# valid primer pairs, scores them, and writes a ranked TSV + diversity PNG.
#
# Dependencies: numpy, matplotlib, Bio (Biopython).  No R required.
#
# Snakemake interface:
#   snakemake.input[0]                  : aligned FASTA (.aln)
#   snakemake.output["tsv"]             : primer pairs TSV
#   snakemake.output["plot"]            : diversity PNG
#   snakemake.params.primer_len         : int
#   snakemake.params.amplicon_min_len   : int
#   snakemake.params.amplicon_max_len   : int
#   snakemake.params.div_cut            : float
#   snakemake.params.GC_tol             : float
#   snakemake.params.min_allele_freq    : float (optional, default 0.05)
#   snakemake.params.max_degeneracy     : int   (optional, default 16)
#   snakemake.log[0]                    : log file
#
# TSV columns (same as R version):
#   primer_id, fwd, rev, fwd_pos, rev_pos, amplicon_len,
#   fwd_GC, rev_GC, fwd_ndegen, rev_ndegen, fwd_fold, rev_fold, total_fold,
#   pair_diversity, delta_GC, combined_score
# =============================================================================

import csv
import logging
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import snakemake
from Bio import AlignIO

# ---------------------------------------------------------------------------
# IUPAC constants
# ---------------------------------------------------------------------------
_IUPAC_COMP_TABLE = str.maketrans(
    "ACGTRYMKSWHBVDNacgtrymkswhbvdn", "TGCAYRKMSWDVBHNtgcayrkmswdvbhn"
)

_IUPAC_MAP = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "AG": "R",
    "CT": "Y",
    "CG": "S",
    "AT": "W",
    "GT": "K",
    "AC": "M",
    "CGT": "B",
    "AGT": "D",
    "ACT": "H",
    "ACG": "V",
    "ACGT": "N",
}

_TSV_FIELDNAMES = [
    "primer_id",
    "fwd",
    "rev",
    "fwd_pos",
    "rev_pos",
    "amplicon_len",
    "fwd_GC",
    "rev_GC",
    "fwd_ndegen",
    "rev_ndegen",
    "fwd_fold",
    "rev_fold",
    "total_fold",
    "pair_diversity",
    "delta_GC",
    "combined_score",
]

_VALID_BASES = frozenset("ACGT")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _rev_comp(seq: str) -> str:
    """IUPAC-aware reverse complement."""
    return seq.translate(_IUPAC_COMP_TABLE)[::-1]


def _iupac_encode(bases: str) -> str:
    """Sorted unique base letters → IUPAC ambiguity code."""
    key = "".join(sorted(set(bases.upper()) & _VALID_BASES))
    return _IUPAC_MAP.get(key, "N")


def _shannon(col: np.ndarray, n_seqs: int) -> float:
    """Shannon entropy of a single alignment column."""
    _, counts = np.unique(col, return_counts=True)
    probs = counts / n_seqs
    probs = probs[probs > 0]
    return float(-np.sum(probs * np.log(probs)))


def _consensus(dna_matrix: np.ndarray) -> np.ndarray:
    """Most frequent base per column."""
    n = dna_matrix.shape[1]
    cons = np.empty(n, dtype="<U1")
    for i in range(n):
        uniq, cnt = np.unique(dna_matrix[:, i], return_counts=True)
        cons[i] = uniq[np.argmax(cnt)]
    return cons


def _rollmean(values: np.ndarray, k: int) -> np.ndarray:
    """Rolling mean with NaN padding at edges (matches R's zoo::rollmean)."""
    result = np.full(len(values), np.nan)
    if k <= len(values):
        conv = np.convolve(values, np.ones(k) / k, mode="valid")
        start = k // 2
        result[start : start + len(conv)] = conv
    return result


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------
def _compute_position_annotations(
    dna_matrix: np.ndarray, n_seqs: int, min_allele_freq: float
):
    """Return (pos_code, pos_fold) — IUPAC letter + fold for each column."""
    aln_len = dna_matrix.shape[1]
    pos_code = np.empty(aln_len, dtype="<U1")
    pos_fold = np.ones(aln_len, dtype=int)

    for i in range(aln_len):
        col = dna_matrix[:, i]
        unique, counts = np.unique(col, return_counts=True)

        # drop gaps and Ns
        mask = ~np.isin(unique, ["-", "n"])
        unique = unique[mask]
        counts = counts[mask]

        if len(unique) == 0:  # all-gap column
            pos_code[i] = "N"
            pos_fold[i] = 1
            continue

        freqs = counts / n_seqs
        keep = unique[freqs >= min_allele_freq]
        keep = np.array([b for b in keep if b.upper() in _VALID_BASES])

        if len(keep) == 0:  # nothing passes freq threshold
            keep = np.array([unique[np.argmax(counts)]])

        pos_code[i] = _iupac_encode("".join(keep))
        pos_fold[i] = len(keep)

    return pos_code, pos_fold


def _build_kmers(consensus, pos_code, pos_fold, divs, primer_len):
    """Slide a window of length *primer_len* across the alignment."""
    n_kmers = len(divs) - primer_len
    kmers = []
    for j in range(n_kmers):
        idx = slice(j, j + primer_len)
        gc_count = int(np.sum(np.isin(consensus[idx], ["g", "c"])))
        kmers.append(
            {
                "pos": j,
                "degen": "".join(pos_code[idx]),
                "n_degen": int(np.sum(pos_fold[idx] > 1)),
                "fold": int(np.prod(pos_fold[idx])),
                "divs": float(np.sum(divs[idx])),
                "GC": gc_count / primer_len,
            }
        )
    return kmers


def _evaluate_pairs(candidates, amplicon_min_len, amplicon_max_len, GC_tol):
    """Return a sorted list of primer-pair dicts (best first)."""
    results = []
    nc = len(candidates)
    for i in range(nc - 1):
        ci = candidates[i]
        for j in range(i + 1, nc):
            cj = candidates[j]
            amp_len = cj["pos"] - ci["pos"]
            if not (amplicon_min_len <= amp_len <= amplicon_max_len):
                continue
            delta_gc = abs(ci["GC"] - cj["GC"])
            if delta_gc >= GC_tol:
                continue

            pair_div = ci["divs"] + cj["divs"]
            total_fold = ci["fold"] * cj["fold"]
            score = 1.0 / (abs(pair_div) + 10.0 * delta_gc**2 + 0.01)

            results.append(
                {
                    "fwd": ci["degen"],
                    "rev": _rev_comp(cj["degen"]),
                    "fwd_pos": ci["pos"],
                    "rev_pos": cj["pos"],
                    "amplicon_len": amp_len,
                    "fwd_GC": round(ci["GC"], 4),
                    "rev_GC": round(cj["GC"], 4),
                    "fwd_ndegen": ci["n_degen"],
                    "rev_ndegen": cj["n_degen"],
                    "fwd_fold": ci["fold"],
                    "rev_fold": cj["fold"],
                    "total_fold": total_fold,
                    "pair_diversity": round(pair_div, 4),
                    "delta_GC": round(delta_gc, 4),
                    "combined_score": round(score, 6),
                }
            )

    results.sort(key=lambda x: (-x["combined_score"], x["total_fold"]))
    for idx, row in enumerate(results, 1):
        row["primer_id"] = f"primer_pair_{idx}"
    return results


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------
def _write_tsv(rows, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=_TSV_FIELDNAMES, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def _write_empty_tsv(path):
    _write_tsv([], path)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def _plot_diversity(
    divs, roll_means, roll_k, results, top_n, primer_len, aln_file, out_plot, log
):
    fig, ax = plt.subplots(figsize=(12, 6))

    ax.scatter(
        np.arange(len(divs)),
        divs,
        s=4,
        alpha=0.5,
        c="black",
        edgecolors="none",
        label="Per-position entropy",
    )

    valid = ~np.isnan(roll_means)
    ax.plot(
        np.arange(len(divs))[valid],
        roll_means[valid],
        color="#2ca25f",
        lw=1.5,
        label=f"Rolling mean (k={roll_k})",
    )

    ymax = float(max(divs)) * 1.05
    for k in range(top_n):
        ax.axvspan(
            results[k]["fwd_pos"],
            results[k]["fwd_pos"] + primer_len,
            alpha=0.20,
            color="#e34a33",
        )
        ax.axvspan(
            results[k]["rev_pos"],
            results[k]["rev_pos"] + primer_len,
            alpha=0.20,
            color="#e34a33",
        )

    from matplotlib.patches import Patch

    handles, _ = ax.get_legend_handles_labels()
    handles.append(
        Patch(facecolor="#e34a33", alpha=0.20, label=f"Top {top_n} primer sites")
    )
    ax.legend(handles=handles, loc="upper right")

    ax.set_title(f"Sequence diversity — {os.path.basename(aln_file)}")
    ax.set_xlabel("Alignment position (bp)")
    ax.set_ylabel("Shannon entropy")
    ax.set_ylim(-0.05, ymax)

    fig.tight_layout()
    fig.savefig(out_plot, dpi=150)
    plt.close(fig)

    log.info("Wrote diversity plot to %s", out_plot)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    # --- Snakemake interface -----------------------------------------------
    aln_file = snakemake.input[0]
    out_tsv = snakemake.output["tsv"]
    out_plot = snakemake.output["plot"]

    log_path = snakemake.log[0]
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    log = logging.getLogger()

    p = snakemake.params
    primer_len = int(p["primer_len"])
    amplicon_min_len = int(p["amplicon_min_len"])
    amplicon_max_len = int(p["amplicon_max_len"])
    div_cut = float(p["div_cut"])
    GC_tol = float(p["GC_tol"])
    min_allele_freq = float(p.get("min_allele_freq", 0.05))
    max_degeneracy = int(p.get("max_degeneracy", 16))

    log.info("Parameters:")
    for k, v in [
        ("aln_file", aln_file),
        ("primer_len", primer_len),
        ("amplicon_min_len", amplicon_min_len),
        ("amplicon_max_len", amplicon_max_len),
        ("div_cut", div_cut),
        ("GC_tol", GC_tol),
        ("min_allele_freq", min_allele_freq),
        ("max_degeneracy", max_degeneracy),
    ]:
        log.info("  %-18s = %s", k, v)

    # --- 1. Load alignment ------------------------------------------------
    aln = AlignIO.read(aln_file, "fasta")
    dna_matrix = np.array([list(str(rec.seq).lower()) for rec in aln])
    n_seqs, aln_len = dna_matrix.shape

    log.info("Loaded %d sequences, alignment length %d bp", n_seqs, aln_len)
    if n_seqs < 2:
        log.error("Need ≥ 2 sequences; only %d found.", n_seqs)
        sys.exit(1)

    # --- 2. Per-position Shannon entropy ----------------------------------
    consensus = _consensus(dna_matrix)
    divs = np.array([_shannon(dna_matrix[:, i], n_seqs) for i in range(aln_len)])
    log.info("Mean per-position entropy: %.4f", np.mean(divs))

    # --- 3. Per-position IUPAC degenerate code + fold ---------------------
    pos_code, pos_fold = _compute_position_annotations(
        dna_matrix, n_seqs, min_allele_freq
    )
    log.info(
        "Degenerate columns (fold > 1): %d / %d", int(np.sum(pos_fold > 1)), aln_len
    )

    # --- 4. Rolling mean for diversity plot -------------------------------
    roll_k = min(10, aln_len)
    roll_means = _rollmean(divs, roll_k)

    # --- 5. Build kmer table ----------------------------------------------
    if aln_len < primer_len:
        log.error("Alignment (%d bp) < primer_len (%d bp).", aln_len, primer_len)
        sys.exit(1)

    kmers = _build_kmers(consensus, pos_code, pos_fold, divs, primer_len)

    # --- 6. Filter candidates --------------------------------------------
    candidates = [
        k for k in kmers if k["divs"] <= div_cut and k["fold"] <= max_degeneracy
    ]

    if not candidates:
        log.info(
            "No candidates pass div_cut=%.2f + max_degeneracy=%d. Empty TSV.",
            div_cut,
            max_degeneracy,
        )
        _write_empty_tsv(out_tsv)
        return

    log.info("%d candidate kmers pass filters", len(candidates))

    # --- 7. Evaluate pairs ------------------------------------------------
    results = _evaluate_pairs(candidates, amplicon_min_len, amplicon_max_len, GC_tol)

    if not results:
        log.info("No valid primer pairs found. Empty TSV.")
        _write_empty_tsv(out_tsv)
        return

    # --- 8. Write TSV -----------------------------------------------------
    _write_tsv(results, out_tsv)
    log.info("Wrote %d primer pairs → %s", len(results), out_tsv)

    # --- 9. Diversity plot ------------------------------------------------
    _plot_diversity(
        divs,
        roll_means,
        roll_k,
        results,
        min(5, len(results)),
        primer_len,
        aln_file,
        out_plot,
        log,
    )


if __name__ == "__main__":
    main()
