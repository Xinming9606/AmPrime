# =============================================================================
# extract_gene.smk
#
# Per-gene rule: extracts target gene sequences from all downloaded CDS and
# RNA FASTA files, producing a single merged FASTA per gene.
#
# Depends on: rule download_genomes. Genome directories and config are passed
# to the CLI as plain paths.
#
# Output : results/{genus}/extracted/{gene}.fasta   [temp]
# =============================================================================

rule extract_gene:
    input:
        cds = str(GENOMES_CDS),
        rna = str(GENOMES_RNA)
    output:
        temp(str(EXTRACTED / "{gene}.fasta"))
    params:
        gene = lambda wc: wc.gene,
        config_file = "config/config.yaml",
        cds_dir = str(GENOMES_CDS),
        rna_dir = str(GENOMES_RNA)
    log:
        str(RESULTS / "logs" / "extract_gene" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "extract_gene" / "{gene}.txt")
    shell:
        """
        python workflow/scripts/extract_gene.py \
            --cds-dir {params.cds_dir:q} \
            --rna-dir {params.rna_dir:q} \
            --out-fasta {output:q} \
            --gene {params.gene:q} \
            --config {params.config_file:q} \
            --log {log:q}
        """
