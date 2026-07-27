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
        config_file = "config/config.yaml"
    log:
        str(RESULTS / "logs" / "check_primers" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "check_primers" / "{gene}.txt")
    shell:
        """
        python workflow/scripts/check_primers.py \
            --in-tsv {input:q} \
            --out-tsv {output:q} \
            --config {params.config_file:q} \
            --log {log:q}
        """
