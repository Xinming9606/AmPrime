#!/usr/bin/env python3
# =============================================================================
# check_primers.py
#
# Primer quality filter: computes hairpin, homodimer, heterodimer, and
# 3'-end stability for every primer pair, then drops pairs that fail any
# configurable threshold.  Writes a filtered TSV with extra quality columns.
#
# Hairpin / dimer ΔG via primer3-py (Santalucia 1998 NN parameters).
# 3'-end ΔG computed directly from the same NN parameters on the terminal
# 5 bases.
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

import primer3
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


def _complement(base: str) -> str:
    return {"A": "T", "T": "A", "C": "G", "G": "C"}.get(base, base)


def calc_3end_dg(seq: str, n_bases: int = 5) -> float:
    """ΔG of the 3'-most *n_bases* of *seq* (kcal/mol, SantaLucia 1998)."""
    end = seq[-n_bases:].upper()
    if len(end) < 2:
        return 0.0

    dg = INIT_GC if end[0] in "GC" else INIT_AT
    for i in range(len(end) - 1):
        dg += NN_DG.get(end[i : i + 2], 0)

    # symmetry correction
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

    hp_f = primer3.calc_hairpin(fwd)
    hp_r = primer3.calc_hairpin(rev)
    hd_f = primer3.calc_homodimer(fwd)
    hd_r = primer3.calc_homodimer(rev)
    het = primer3.calc_heterodimer(fwd, rev)
    dg3_f = calc_3end_dg(fwd)
    dg3_r = calc_3end_dg(rev)

    metrics = {
        "hairpin_fwd_dg": round(hp_f.dg, 2),
        "hairpin_rev_dg": round(hp_r.dg, 2),
        "homodimer_fwd_dg": round(hd_f.dg, 2),
        "homodimer_rev_dg": round(hd_r.dg, 2),
        "heterodimer_dg": round(het.dg, 2),
        "end3_fwd_dg": dg3_f,
        "end3_rev_dg": dg3_r,
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
