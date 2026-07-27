# =============================================================================
# reports.smk
#
# gene_report : per-gene HTML report (primers + PCR + diversity)
# =============================================================================


rule gene_report:
    input:
        primers   = str(PRIMERS / "{gene}_primers.tsv"),
        amplicons = str(PRIMERS / "{gene}_amplicons.tsv"),
        diversity = str(PRIMERS / "{gene}_diversity.png")
    output:
        str(REPORTS / "{gene}_report.html")
    params:
        gene  = lambda wc: wc.gene,
        genus = config["genus"]
    log:
        str(RESULTS / "logs" / "gene_report" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "gene_report" / "{gene}.txt")
    shell:
        """
        python workflow/scripts/gene_report.py \
            --gene {params.gene:q} \
            --genus {params.genus:q} \
            --primers-tsv {input.primers:q} \
            --amplicons-tsv {input.amplicons:q} \
            --diversity-png {input.diversity:q} \
            --out-html {output:q} \
            --log {log:q}
        """
