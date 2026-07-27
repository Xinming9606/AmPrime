#!/usr/bin/env python3
# =============================================================================
# check_primers.py
#
# Primer quality filter: computes hairpin, homodimer, heterodimer, and
# 3'-end stability for every primer pair, then drops pairs that fail any
# configurable threshold.  Writes a filtered TSV with extra quality columns.
#
# All ΔG calculations use SantaLucia 1998 nearest-neighbour parameters in
# pure Python — no primer3 C library required.  Works on any OS.
#
# Snakemake interface:
#   snakemake.input[0]             : {gene}_primers_raw.tsv  (from design_primers)
#   snakemake.output[0]            : {gene}_primers.tsv       (filtered + columns)
#   snakemake.params.max_hairpin_dg       : float or None
#   snakemake.params.max_homodimer_dg     : float or None
#   snakemake.params.max_heterodimer_dg   : float or None
#   snakemake.params.max_3end_dg          : float or None
#   snakemake.log[0]               : log file
# =============================================================================

import csv
import logging
import os

import snakemake

# ---------------------------------------------------------------------------
# SantaLucia 1998 unified nearest-neighbour ΔG (kcal/mol, 37 °C, 1 M NaCl)
# ---------------------------------------------------------------------------
NN_DG = {
    "AA": -1.00,
    "TT": -1.00,
    "AT": -0.88,
    "TA": -0.58,
    "CA": -1.45,
    "TG": -1.45,
    "GT": -1.44,
    "AC": -1.44,
    "CT": -1.28,
    "AG": -1.28,
    "GA": -1.30,
    "TC": -1.30,
    "CG": -2.17,
    "GC": -2.24,
    "GG": -1.84,
    "CC": -1.84,
}

SYMMETRY_DG = 0.43  # self-complementary duplex penalty
INIT_AT = 2.30  # initiation — terminal A·T
INIT_GC = 2.55  # initiation — terminal G·C
MIN_STEM = 2  # minimum stem length for hairpin
MIN_LOOP = 3  # minimum loop size for hairpin


# ---------------------------------------------------------------------------
# Building blocks
# ---------------------------------------------------------------------------
def _complement(base: str) -> str:
    return {"A": "T", "T": "A", "C": "G", "G": "C"}.get(base, base)


def _reverse_complement(seq: str) -> str:
    return "".join(_complement(b) for b in reversed(seq))


def _duplex_dg(seq_5p: str) -> float:
    """ΔG of a duplex scored from one strand (5'→3'), SantaLucia 1998.

    seq_5p is the sequence of one strand in 5'→3' orientation.  The
    complementary strand is assumed to pair perfectly.
    """
    n = len(seq_5p)
    if n < 2:
        return 0.0
    dg = INIT_GC if seq_5p[0] in "GC" else INIT_AT
    for i in range(n - 1):
        dg += NN_DG.get(seq_5p[i : i + 2], 0)
    if seq_5p == _reverse_complement(seq_5p):
        dg += SYMMETRY_DG
    return dg


# ---------------------------------------------------------------------------
# Dimer ΔG (homo- and hetero-)
# ---------------------------------------------------------------------------
def _align_dg(seq_a: str, seq_b_rc: str) -> float:
    """Best (most negative) ΔG of antiparallel alignment of seq_a vs seq_b_rc.

    seq_b_rc is the reverse complement of some other sequence in 5'→3'.
    We slide the two 5'→3' strands across each other and score every
    overlapping region.
    """
    best = 0.0
    la, lb = len(seq_a), len(seq_b_rc)

    for offset in range(-(lb - 1), la):
        if offset >= 0:
            a0, _, n = offset, 0, min(la - offset, lb)
        else:
            a0, _, n = 0, -offset, min(la, lb + offset)

        if n < MIN_STEM:
            continue

        dg = _duplex_dg(seq_a[a0 : a0 + n])
        best = min(best, dg)

    return round(best, 2)


def calc_homodimer_dg(seq: str) -> float:
    """Minimum homodimer ΔG (self-dimer)."""
    return _align_dg(seq, _reverse_complement(seq))


def calc_heterodimer_dg(seq_a: str, seq_b: str) -> float:
    """Minimum heterodimer ΔG (cross-dimer)."""
    return _align_dg(seq_a, _reverse_complement(seq_b))


# ---------------------------------------------------------------------------
# Hairpin ΔG
# ---------------------------------------------------------------------------
def calc_hairpin_dg(seq: str) -> float:
    """Minimum hairpin ΔG — find the most stable stem within the sequence.

    Enumerates all possible stem-start (i), stem-end-start (j) and stem
    length (k) combinations with a minimum loop of MIN_LOOP bases between the
    two stem halves.  Because primer oligos are short (15–35 nt) the
    triple-loop is harmless.
    """
    n = len(seq)
    best = 0.0

    for stem_len in range(MIN_STEM, n // 2 + 1):
        max_i = n - 2 * stem_len - MIN_LOOP
        for i in range(max_i + 1):
            j_min = i + stem_len + MIN_LOOP
            for _ in range(j_min, n - stem_len + 1):
                stem5 = seq[i : i + stem_len]
                dg = _duplex_dg(stem5)
                best = min(best, dg)

    return round(best, 2)


# ---------------------------------------------------------------------------
# 3'-end stability
# ---------------------------------------------------------------------------
def calc_3end_dg(seq: str, n_bases: int = 5) -> float:
    """ΔG of the 3'-most *n_bases* of *seq* (kcal/mol, SantaLucia 1998)."""
    end = seq[-n_bases:].upper()
    if len(end) < 2:
        return 0.0

    dg = INIT_GC if end[0] in "GC" else INIT_AT
    for i in range(len(end) - 1):
        dg += NN_DG.get(end[i : i + 2], 0)

    rc = "".join(_complement(b) for b in reversed(end))
    if end == rc:
        dg += SYMMETRY_DG

    return round(dg, 2)


# ---------------------------------------------------------------------------
# Quality-check a single primer pair
# ---------------------------------------------------------------------------
def check_pair(row: dict, thresholds: dict) -> dict:
    fwd = row["fwd"]
    rev = row["rev"]

    metrics = {
        "hairpin_fwd_dg": calc_hairpin_dg(fwd),
        "hairpin_rev_dg": calc_hairpin_dg(rev),
        "homodimer_fwd_dg": calc_homodimer_dg(fwd),
        "homodimer_rev_dg": calc_homodimer_dg(rev),
        "heterodimer_dg": calc_heterodimer_dg(fwd, rev),
        "end3_fwd_dg": calc_3end_dg(fwd),
        "end3_rev_dg": calc_3end_dg(rev),
    }

    # thresholds: more negative ΔG = stronger secondary structure = worse
    checks = [
        ("hairpin_fwd_dg", thresholds.get("max_hairpin_dg")),
        ("hairpin_rev_dg", thresholds.get("max_hairpin_dg")),
        ("homodimer_fwd_dg", thresholds.get("max_homodimer_dg")),
        ("homodimer_rev_dg", thresholds.get("max_homodimer_dg")),
        ("heterodimer_dg", thresholds.get("max_heterodimer_dg")),
        ("end3_fwd_dg", thresholds.get("max_3end_dg")),
        ("end3_rev_dg", thresholds.get("max_3end_dg")),
    ]

    fails = []
    for key, limit in checks:
        if limit is not None and metrics[key] < limit:
            fails.append(key)

    metrics["qc_pass"] = len(fails) == 0
    metrics["qc_fail_reasons"] = ";".join(fails) if fails else ""

    return metrics


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    log_path = snakemake.log[0]
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    log = logging.getLogger()

    in_tsv = snakemake.input[0]
    out_tsv = snakemake.output[0]

    # helper: cast params to float or None
    def p(name):
        v = snakemake.params.get(name)
        try:
            return float(v) if v is not None else None
        except (TypeError, ValueError):
            return None

    thresholds = {
        "max_hairpin_dg": p("max_hairpin_dg"),
        "max_homodimer_dg": p("max_homodimer_dg"),
        "max_heterodimer_dg": p("max_heterodimer_dg"),
        "max_3end_dg": p("max_3end_dg"),
    }

    log.info("Input  : %s", in_tsv)
    log.info("Output : %s", out_tsv)
    log.info("Thresholds: %s", {k: v for k, v in thresholds.items() if v is not None})

    # Read input
    with open(in_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)

    if not rows:
        log.warning("Input TSV is empty — writing empty output.")
        os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
        with open(out_tsv, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            # minimal header when there's no data
            if reader.fieldnames:
                extra = [
                    "hairpin_fwd_dg",
                    "hairpin_rev_dg",
                    "homodimer_fwd_dg",
                    "homodimer_rev_dg",
                    "heterodimer_dg",
                    "end3_fwd_dg",
                    "end3_rev_dg",
                    "qc_pass",
                    "qc_fail_reasons",
                ]
                w.writerow(reader.fieldnames + extra)
        return

    in_fields = reader.fieldnames or []
    qc_fields = [
        "hairpin_fwd_dg",
        "hairpin_rev_dg",
        "homodimer_fwd_dg",
        "homodimer_rev_dg",
        "heterodimer_dg",
        "end3_fwd_dg",
        "end3_rev_dg",
        "qc_pass",
        "qc_fail_reasons",
    ]
    out_fields = in_fields + qc_fields

    checked = []
    for r in rows:
        m = check_pair(r, thresholds)
        checked.append({**r, **m})

    passed = [r for r in checked if r["qc_pass"]]
    n_dropped = len(checked) - len(passed)
    log.info(
        "Checked %d pairs — %d passed, %d dropped.",
        len(checked),
        len(passed),
        n_dropped,
    )

    if not passed:
        log.warning(
            "All %d pairs failed quality filters. Writing empty TSV.", len(checked)
        )
        os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
        with open(out_tsv, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(out_fields)
        return

    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=out_fields, delimiter="\t", extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(passed)

    log.info("Written %d filtered pairs to %s", len(passed), out_tsv)


if __name__ == "__main__":
    main()
