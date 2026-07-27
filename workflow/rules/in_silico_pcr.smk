# =============================================================================
# in_silico_pcr.smk
#
# Per-gene rule: validates the top-ranked QC-passed primer candidates against
# the full genome assemblies using the Python primer scanner. Reports are sorted
# by in silico PCR performance so the HTML report can recommend the best hit.
#
# Inputs : results/{genus}/primers/{gene}_primers.tsv   (from design_primers)
#          results/{genus}/genomes/genomic/             (from download_genomes)
# Output : results/{genus}/primers/{gene}_amplicons.tsv
# =============================================================================

rule in_silico_pcr:
    input:
        primers    = str(PRIMERS / "{gene}_primers.tsv"),
        genome_dir = str(GENOMES_GENOMIC)
    output:
        str(PRIMERS / "{gene}_amplicons.tsv")
    params:
        gene             = lambda wc: wc.gene,
        config_file = "config/config.yaml"
    log:
        str(RESULTS / "logs" / "in_silico_pcr" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "in_silico_pcr" / "{gene}.txt")
    shell:
        """
        python workflow/scripts/in_silico_pcr.py \
            --primers-tsv {input.primers:q} \
            --genome-dir {input.genome_dir:q} \
            --out-tsv {output:q} \
            --gene {params.gene:q} \
            --config {params.config_file:q} \
            --log {log:q}
        """
