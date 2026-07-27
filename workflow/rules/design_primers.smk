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
        primer_len       = config["primer_len"],
        amplicon_min_len = config["amplicon_min_len"],
        amplicon_max_len = config["amplicon_max_len"],
        div_cut          = lambda wc: config.get("div_cut_per_gene", {}).get(wc.gene, config["div_cut"]),
        GC_tol           = config["GC_tol"],
        min_allele_freq  = config.get("min_allele_freq", 0.05),
        max_degeneracy   = config.get("max_degeneracy", 16)
    log:
        str(RESULTS / "logs" / "design_primers" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "design_primers" / "{gene}.txt")
    shell:
        """
        python workflow/scripts/design_primers.py \
            --aln {input.aln:q} \
            --out-tsv {output.tsv:q} \
            --out-plot {output.plot:q} \
            --primer-len {params.primer_len} \
            --amplicon-min-len {params.amplicon_min_len} \
            --amplicon-max-len {params.amplicon_max_len} \
            --div-cut {params.div_cut} \
            --gc-tol {params.GC_tol} \
            --min-allele-freq {params.min_allele_freq} \
            --max-degeneracy {params.max_degeneracy} \
            --log {log:q}
        """
