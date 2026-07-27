# =============================================================================
# align.smk
#
# Per-gene rule: multiple sequence alignment of the clustered centroids. The
# default backend is Python for cross-platform installs; config.yaml can select
# auto/MAFFT/MUSCLE when external aligners are available.
#
# Input  : results/{genus}/extracted/{gene}.centroids.fasta  [temp, from cluster]
# Output : results/{genus}/aligned/{gene}.aln                [temp]
#
# .aln is temp(): consumed only by design_primers; not needed afterwards.
# =============================================================================


rule align:
    input:
        str(EXTRACTED / "{gene}.centroids.fasta")
    output:
        temp(str(ALIGNED / "{gene}.aln"))
    params:
        config_file = "config/config.yaml"
    log:
        str(RESULTS / "logs" / "align" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "align" / "{gene}.txt")
    shell:
        """
        python workflow/scripts/align_fasta.py \
            --input {input:q} \
            --output {output:q} \
            --config {params.config_file:q} \
            --log {log:q}
        """
