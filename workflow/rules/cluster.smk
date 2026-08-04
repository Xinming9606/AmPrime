# =============================================================================
# cluster.smk
#
# Per-gene rule: dereplicates near-identical sequences with a Python greedy
# centroid clusterer at 97% identity, keeping one centroid per cluster.
# This reduces redundancy before alignment.
#
# Input  : results/{genus}/extracted/{gene}.fasta            [temp, from extract_gene]
# Output : results/{genus}/extracted/{gene}.centroids.fasta  [temp]
#
# centroids.fasta is temp(): it feeds align only. in_silico_pcr validates
# against the full genomes (genomes/genomic/), not the centroids.
# =============================================================================

rule cluster:
    input:
        fasta = str(EXTRACTED / "{gene}.fasta"),
        config = CONFIG_FILE
    output:
        temp(str(EXTRACTED / "{gene}.centroids.fasta"))
    params:
        identity = 0.97,
        script = str(SCRIPTS_DIR / "cluster_fasta.py")
    log:
        str(RESULTS / "logs" / "cluster" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "cluster" / "{gene}.txt")
    shell:
        """
        python {params.script:q} \
            --input {input.fasta:q} \
            --output {output:q} \
            --identity {params.identity} \
            --log {log:q}
        """
