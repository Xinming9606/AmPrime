# AmPrime

AmPrime designs and validates amplicon-sequencing primer pairs for bacterial
housekeeping genes using public NCBI genomes.

You give it:

- a bacterial genus, such as `Borrelia`
- one or more target genes, such as `recG`, `clpA`, or an MLST gene set

It returns one self-contained HTML report per gene with the recommended primer
pair, in silico PCR validation, sequence diversity plot, and all candidate
primer pairs.

## When To Use

Use AmPrime when you want a reproducible first-pass primer design workflow for
a bacterial genus, especially when you want primers that should work across
many genomes within that genus.

It is not a full specificity checker yet. The current pipeline checks whether
the top primer pair amplifies genomes inside the target genus, but it does not
test off-target amplification outside the genus.

## How It Works

For each gene independently, AmPrime runs this pipeline:

```mermaid
flowchart LR
    A[Download genomes<br/>NCBI] --> B[Extract gene<br/>CDS + rRNA]
    B --> C[Dereplicate<br/>vsearch 97%]
    C --> D[Align<br/>MUSCLE]
    D --> E[Design primers<br/>entropy scan]
    E --> F[Quality filter<br/>hairpin + dimer]
    F --> G[In silico PCR<br/>seqkit amplicon]
    G --> H[HTML report]
```

The main idea is simple:

1. Download genomic, CDS, and RNA FASTA files for the target genus.
2. Extract the target gene from CDS/RNA annotations.
3. Dereplicate near-identical sequences so redundant strains do not dominate.
4. Align representative sequences.
5. Find conserved primer windows that flank a variable amplicon region.
6. Filter primer pairs for simple secondary-structure risks.
7. Validate the top pair against full genomes with `seqkit amplicon`.
8. Write a browsable HTML report.

Missing genes are handled gracefully. If a gene cannot be found, the pipeline
continues and writes a report showing that no candidates were available.

## Quick Start

```bash
# 1. Clone the repository
git clone git@github.com:Xinming9606/AmPrime.git
cd AmPrime

# 2. Create and activate the environment
micromamba env create -f workflow/envs/environment.yaml
micromamba activate primer-pipeline

# 3. Edit config/config.yaml
# Set your genus, genes, and optional gene aliases.

# 4. Preview the work to be done
snakemake -n

# 5. Run the pipeline
snakemake --cores 4
```

Reports are written to:

```text
results/<genus>/reports/<gene>_report.html
```

Open the HTML file in a browser to inspect the recommended primer pair and
validation summary.

## Configuration

All user-facing settings live in [config/config.yaml](config/config.yaml).

Minimal example:

```yaml
genus: Borrelia

genes:
  - recG
  - clpA
  - uvrA

assembly_level: complete

primer_len: 20
amplicon_min_len: 300
amplicon_max_len: 1000
div_cut: 2.0
GC_tol: 0.1

pcr_mismatch: 3
```

Useful options:

| Setting | Meaning |
| --- | --- |
| `genus` | Bacterial genus name recognized by NCBI. |
| `genes` | One or more gene names. Each gene is processed independently. |
| `assembly_level` | NCBI assembly level: `complete`, `chromosome`, `scaffold`, or `contig`. |
| `primer_len` | Primer length in bp. |
| `amplicon_min_len`, `amplicon_max_len` | Target amplicon size range. |
| `div_cut` | Maximum Shannon entropy allowed for conserved primer windows. Raise it if no primers are found. |
| `GC_tol` | Maximum GC fraction difference between forward and reverse primers. |
| `pcr_mismatch` | Mismatches allowed by `seqkit amplicon`. |

If a gene is annotated under different names across genomes, add aliases:

```yaml
gene_aliases:
  tuf:
    - tsf
  16S:
    - "16S ribosomal RNA"
```

You can also override the diversity cutoff for specific genes:

```yaml
div_cut_per_gene:
  16S: 3.0
  rpoB: 1.5
```

## Output Files

For each gene, outputs are written under `results/<genus>/`.

| File | Contents |
| --- | --- |
| `reports/<gene>_report.html` | Main deliverable: recommendation, PCR validation, plot, and candidate table. |
| `primers/<gene>_primers.tsv` | Filtered candidate primer pairs ranked by score. |
| `primers/<gene>_amplicons.tsv` | In silico PCR result for the top primer pair. |
| `primers/<gene>_diversity.png` | Per-position Shannon entropy plot with top primer sites marked. |
| `logs/...` | Per-step logs for debugging. |
| `benchmarks/...` | Snakemake benchmark files. |

## Project Layout

```text
AmPrime/
|-- Snakefile
|-- config/
|   `-- config.yaml
|-- workflow/
|   |-- Snakefile
|   |-- rules/
|   |   |-- download_genomes.smk
|   |   |-- extract_gene.smk
|   |   |-- cluster.smk
|   |   |-- align.smk
|   |   |-- design_primers.smk
|   |   |-- check_primers.smk
|   |   |-- in_silico_pcr.smk
|   |   `-- reports.smk
|   |-- scripts/
|   |   |-- download_genomes.py
|   |   |-- extract_gene.py
|   |   |-- cluster_fasta.py
|   |   |-- align_fasta.py
|   |   |-- design_primers.py
|   |   |-- check_primers.py
|   |   |-- in_silico_pcr.py
|   |   |-- gene_report.py
|   |   `-- gene_report.html
|   `-- envs/
|       `-- environment.yaml
`-- results/
```

## Design Notes

Snakemake is intentionally kept as a thin scheduler. It manages dependencies,
parallel execution, logs, benchmarks, and resumability. The actual work is done
by standalone Python command-line tools in `workflow/scripts/`.

This keeps each step easy to test and debug. For example:

```bash
python workflow/scripts/extract_gene.py --help
python workflow/scripts/design_primers.py --help
python workflow/scripts/gene_report.py --help
```

## Troubleshooting

If the dry run fails, check that you are running from the repository root:

```bash
snakemake -n
```

If no primers are found, try one or more of the following:

- raise `div_cut`
- relax `GC_tol`
- add gene aliases under `gene_aliases`
- use a less restrictive `assembly_level`
- inspect the per-gene logs under `results/<genus>/logs/`

If genome download fails, confirm that the genus name is recognized by NCBI and
that your internet connection is available.

If `MUSCLE` is slow, start with a small genus or a stricter assembly level such
as `complete`.

## Requirements

You need a conda-compatible package manager such as `conda`, `mamba`, or
`micromamba`, plus an internet connection for the genome download step.

All workflow dependencies are pinned in:

```text
workflow/envs/environment.yaml
```

The environment includes Snakemake, Python, Biopython, NumPy, Matplotlib,
`ncbi-genome-download`, `vsearch`, `MUSCLE`, `seqkit`, and Python `markdown`.

## Limitations

- Off-target specificity outside the target genus is not checked yet.
- Primer windows are derived from the alignment consensus.
- Degenerate-base handling is conservative.
- Very large genera can take a long time to download and align.

## License

MIT. See [LICENSE](LICENSE).
