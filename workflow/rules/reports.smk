# =============================================================================
# reports.smk
#
# Per-gene HTML reports plus a cross-gene comparison report.
# =============================================================================


rule gene_report:
    input:
        primers = str(PRIMERS / "{gene}_primers.tsv"),
        amplicons = str(PRIMERS / "{gene}_amplicons.tsv"),
        species_summary = str(PRIMERS / "{gene}_species_summary.tsv"),
        diversity = str(PRIMERS / "{gene}_diversity.png"),
        alignment_meta = str(ALIGNED / "{gene}.alignment.tsv"),
        config = CONFIG_FILE
    output:
        str(REPORTS / "{gene}_report.html")
    params:
        gene = lambda wc: wc.gene,
        config_file = CONFIG_FILE,
        download_manifest = str(RESULTS / "genomes" / "download_manifest.tsv"),
        script = str(SCRIPTS_DIR / "gene_report.py")
    log:
        str(RESULTS / "logs" / "gene_report" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "gene_report" / "{gene}.txt")
    shell:
        """
        python {params.script:q} \
            --gene {params.gene:q} \
            --config {params.config_file:q} \
            --primers-tsv {input.primers:q} \
            --amplicons-tsv {input.amplicons:q} \
            --species-summary {input.species_summary:q} \
            --diversity-png {input.diversity:q} \
            --alignment-meta {input.alignment_meta:q} \
            --download-manifest {params.download_manifest:q} \
            --out-html {output:q} \
            --log {log:q}
        """


rule cross_gene_report:
    input:
        summaries = expand(str(PRIMERS / "{gene}_species_summary.tsv"), gene=GENES),
        config = CONFIG_FILE
    output:
        str(REPORTS / "cross_gene_report.html")
    params:
        script = str(SCRIPTS_DIR / "cross_gene_report.py")
    log:
        str(RESULTS / "logs" / "cross_gene_report.log")
    benchmark:
        str(RESULTS / "benchmarks" / "cross_gene_report.txt")
    shell:
        """
        python {params.script:q} \
            --summary-tsv {input.summaries:q} \
            --out-html {output:q}
        """
