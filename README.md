# AmPrime

> Design and validate amplicon-sequencing primers for **any bacterial genus**
> and **any set of housekeeping genes** — straight from public NCBI genomes.

A reproducible primer-design pipeline that turns a genus name and a list of
genes into ranked, in-silico-validated primer pairs, each with a self-contained
HTML report. Snakemake is used as a thin scheduler; the per-step Python scripts
are ordinary command-line tools.

## How it works

Give it a genus (e.g. `Borrelia`) and one or more genes (e.g. the eight-gene
MLST scheme). For **each gene independently**, the pipeline:

```mermaid
flowchart LR
    A[Download genomes<br/>NCBI] --> B[Extract gene<br/>CDS + rRNA]
    B --> C[Dereplicate<br/>vsearch 97%]
    C --> D[Align<br/>MUSCLE]
    D --> E[Design primers<br/>Shannon entropy scan]
    E --> F[In silico PCR<br/>seqkit amplicon]
    F --> G[HTML report]
```

1. **Download** every assembly for the genus (`ncbi-genome-download`), keeping
   the genomic, CDS, and rRNA FASTA for each.
2. **Extract** the target gene from the CDS / rRNA FASTA of every genome.
   Missing genes are warned about and skipped, never fatal.
3. **Dereplicate** near-identical sequences at 97% identity (`vsearch`) so the
   alignment is not dominated by redundant strains.
4. **Align** the representatives (`MUSCLE`).
5. **Design primers** by scanning the alignment for conserved windows (low
   Shannon entropy) that flank a variable region of the right amplicon size,
   scoring each pair on conservation and GC balance.
6. **Validate** the top pair by in silico PCR against the *full* genome set,
   reporting how many genomes amplify and the product size.
7. **Report** everything as a single browsable HTML file per gene.

Genes run in parallel — point it at eight genes and it builds eight reports.

## Quick start

```bash
# 1. Clone
git clone git@github.com:Xinming9606/AmPrime.git
cd AmPrime

# 2. Create the environment (one env, all tools)
micromamba env create -f workflow/envs/environment.yaml
micromamba activate primer-pipeline

# 3. Edit config/config.yaml — set your genus and genes

# 4. Check the plan, then run
snakemake -n            # dry-run: shows the DAG without executing
snakemake --cores 4     # real run
```

Outputs land in `results/<genus>/reports/<gene>_report.html`.

## Configuration

Everything is controlled from `config/config.yaml`:

```yaml
genus: Borrelia          # any NCBI bacterial genus
genes:                   # one or more; each runs independently
  - clpA
  - recG
  - uvrA

assembly_level: complete # complete | chromosome | scaffold | contig

primer_len: 20           # primer length (bp)
amplicon_min_len: 300    # target amplicon size range (bp)
amplicon_max_len: 1000
div_cut: 2.0             # max Shannon entropy for a conserved primer window
                         # (raise if no primers are found)
GC_tol: 0.1              # max GC% difference within a primer pair

pcr_mismatch: 3          # mismatches allowed in in silico PCR
```

If a gene's annotation name varies across genomes, add aliases:

```yaml
gene_aliases:
  tuf:
    - tsf
  16S:
    - "16S ribosomal RNA"
```

## Output

For each gene, under `results/<genus>/`:

| File                           | What it contains                                                                |
| ------------------------------ | ------------------------------------------------------------------------------- |
| `reports/<gene>_report.html`   | The deliverable: top pair, PCR validation, diversity plot, full candidate table |
| `primers/<gene>_primers.tsv`   | All candidate primer pairs, ranked by score                                     |
| `primers/<gene>_amplicons.tsv` | In silico PCR result for the top pair                                           |
| `primers/<gene>_diversity.png` | Per-position entropy with primer sites marked                                   |

## Repository layout

```text
AmPrime/
|-- config/
|   `-- config.yaml              # the only file you edit
|-- workflow/
|   |-- Snakefile                # orchestrator: config, paths, includes
|   |-- rules/                   # one .smk module per step
|   |   |-- download_genomes.smk
|   |   |-- extract_gene.smk
|   |   |-- cluster.smk
|   |   |-- align.smk
|   |   |-- design_primers.smk
|   |   |-- in_silico_pcr.smk
|   |   `-- reports.smk
|   |-- scripts/                 # standalone CLIs called by the rules
|   |   |-- extract_gene.py
|   |   |-- design_primers.py
|   |   |-- check_primers.py
|   |   |-- in_silico_pcr.py
|   |   |-- gene_report.py
|   |   `-- gene_report.html
|   `-- envs/
|       `-- environment.yaml     # all dependencies, one environment
`-- results/                     # generated output (git-ignored)
```

## Requirements

A conda-compatible package manager (`conda`, `mamba`, or `micromamba`) and an
internet connection for the genome download step. Everything else - Snakemake,
`ncbi-genome-download`, `vsearch`, `MUSCLE`, `seqkit`, and Python reporting
packages - is pinned in `workflow/envs/environment.yaml`.

## Notes & limitations

- **Specificity against off-target taxa is not yet checked.** The pipeline
  confirms primers amplify within the target genus; it does not test whether
  they also amplify outside it.
- Primer windows are taken from the consensus; degenerate-base handling is
  conservative.
- `MUSCLE` can be slow on genera with thousands of assemblies. Start with a
  small genus to validate your settings.

## License

MIT — see [LICENSE](LICENSE).
