# =============================================================================
# check_primers.smk
#
# Per-gene rule: computes hairpin, dimer, and 3'-end stability metrics for
# every primer pair, filters out pairs that fail configurable thresholds,
# and writes the result as the final primer TSV consumed by downstream rules.
#
# Input  : results/{genus}/primers/{gene}_primers_raw.tsv   [temp, from design_primers]
# Output : results/{genus}/primers/{gene}_primers.tsv       [final]
# =============================================================================


rule check_primers:
    input:
        str(PRIMERS / "{gene}_primers_raw.tsv")
    output:
        str(PRIMERS / "{gene}_primers.tsv")
    params:
        max_hairpin_dg    = config.get("max_hairpin_dg"),
        max_homodimer_dg  = config.get("max_homodimer_dg"),
        max_heterodimer_dg = config.get("max_heterodimer_dg"),
        max_3end_dg        = config.get("max_3end_dg"),
    log:
        str(RESULTS / "logs" / "check_primers" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "check_primers" / "{gene}.txt")
    shell:
        """
        python workflow/scripts/check_primers.py \
            --in-tsv {input:q} \
            --out-tsv {output:q} \
            --max-hairpin-dg {params.max_hairpin_dg} \
            --max-homodimer-dg {params.max_homodimer_dg} \
            --max-heterodimer-dg {params.max_heterodimer_dg} \
            --max-3end-dg {params.max_3end_dg} \
            --log {log:q}
        """
