# =============================================================================
# download_genomes.smk
#
# Downloads three FASTA types for the target genus from NCBI, each into its
# own subdirectory:
#   genomic   - full genome assemblies, used by in_silico_pcr (KEPT)
#   cds       - CDS FASTA, used to extract protein-coding genes
#   rna       - RNA FASTA, used to extract rRNA genes
#
# Uses ncbi-genome-download with --genera (no taxid resolution needed).
#
# Outputs:
#   results/{genus}/genomes/genomic/   - *_genomic.fna           [KEPT]
#   results/{genus}/genomes/cds/       - *_cds_from_genomic.fna
#   results/{genus}/genomes/rna/       - *_rna_from_genomic.fna
#   results/{genus}/genomes/download_manifest.tsv
# =============================================================================


rule download_genomes:
    output:
        genomic = directory(str(GENOMES_GENOMIC)),
        cds     = directory(str(GENOMES_CDS)),
        rna     = directory(str(GENOMES_RNA)),
        manifest = str(RESULTS / "genomes" / "download_manifest.tsv")
    params:
        config_file = "config/config.yaml",
        script = str(SCRIPTS_DIR / "download_genomes.py")
    log:
        str(RESULTS / "logs" / "download_genomes.log")
    benchmark:
        str(RESULTS / "benchmarks" / "download_genomes.txt")
    shell:
        """
        python {params.script:q} \
            --config {params.config_file:q} \
            --genomic-dir {output.genomic:q} \
            --cds-dir {output.cds:q} \
            --rna-dir {output.rna:q} \
            --manifest {output.manifest:q} \
            --log {log:q}
        """
