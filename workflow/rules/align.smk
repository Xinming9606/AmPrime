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
        mkdir -p $(dirname {output})

        n_in=$(grep -c "^>" {input} 2>/dev/null || true)
        n_in=${{n_in:-0}}
        if [ "$n_in" -lt 2 ]; then
            cp {input} {output}
            echo "Only $n_in centroid sequence(s); skipped MUSCLE alignment" >> {log}
            exit 0
        fi

        muscle -align {input} -output {output} 2>> {log}

        n=$(grep -c "^>" {output} 2>/dev/null || true)
        n=${{n:-0}}
        echo "Aligned $n sequences" >> {log}
        """
