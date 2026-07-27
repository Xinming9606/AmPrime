# =============================================================================
# align.smk
#
# Per-gene rule: multiple sequence alignment of the clustered centroids
# using MUSCLE. The alignment feeds design_primers.py for consensus and
# Shannon-entropy calculation.
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
    log:
        str(RESULTS / "logs" / "align" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "align" / "{gene}.txt")
    shell:
        """
        python workflow/scripts/align_fasta.py \
            --input {input:q} \
            --output {output:q} \
            --log {log:q}
        """
