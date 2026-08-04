# =============================================================================
# design_primers.smk
#
# Per-gene rule: runs design_primers.py on the aligned FASTA to produce a
# ranked primer-pair TSV and a diversity PNG with primer sites overlaid.
#
# Inputs  : results/{genus}/aligned/{gene}.aln
# Outputs : results/{genus}/primers/{gene}_primers_raw.tsv  [temp, feeds check_primers]
#           results/{genus}/primers/{gene}_diversity.png
# =============================================================================


rule design_primers:
    input:
        aln = str(RESULTS / "aligned" / "{gene}.aln")
    output:
        tsv  = temp(str(PRIMERS / "{gene}_primers_raw.tsv")),
        plot = str(PRIMERS / "{gene}_diversity.png")
    params:
        gene = lambda wc: wc.gene,
        config_file = "config/config.yaml",
        script = str(SCRIPTS_DIR / "design_primers.py")
    log:
        str(RESULTS / "logs" / "design_primers" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "design_primers" / "{gene}.txt")
    shell:
        """
        python {params.script:q} \
            --aln {input.aln:q} \
            --out-tsv {output.tsv:q} \
            --out-plot {output.plot:q} \
            --gene {params.gene:q} \
            --config {params.config_file:q} \
            --log {log:q}
        """
