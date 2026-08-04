# =============================================================================
# align.smk
#
# Per-gene rule: multiple sequence alignment of the clustered centroids. The
# default backend is Python for cross-platform installs; config.yaml can select
# auto/MAFFT/MUSCLE when external aligners are available.
#
# Input  : results/{genus}/extracted/{gene}.centroids.fasta  [temp, from cluster]
# Output : results/{genus}/aligned/{gene}.aln                [temp]
#          results/{genus}/aligned/{gene}.alignment.tsv      [backend metadata]
#
# .aln is temp(): consumed only by primers_design; not needed afterwards.
# =============================================================================


rule align:
    threads: 1
    resources:
        mem_mb=4096
    input:
        centroids = str(EXTRACTED / "{gene}.centroids.fasta"),
        config = CONFIG_FILE
    output:
        aln = temp(str(ALIGNED / "{gene}.aln")),
        metadata = str(ALIGNED / "{gene}.alignment.tsv")
    params:
        config_file = CONFIG_FILE,
        script = str(SCRIPTS_DIR / "fasta_align.py")
    log:
        str(RESULTS / "logs" / "align" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "align" / "{gene}.txt")
    shell:
        """
        python {params.script:q} \
            --input {input.centroids:q} \
            --output {output.aln:q} \
            --metadata {output.metadata:q} \
            --config {params.config_file:q} \
            --log {log:q}
        """
