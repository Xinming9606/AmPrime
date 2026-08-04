# Functional Test

Use the AmPrime API to run the Borrelia dataset. CI downloads the archive from
NCBI through the project downloader before running the same test.

## Step 1. Check The Dataset

For the full local test, the archive should exist:

```text
data/borrelia-genomes.tar.gz
```

It should contain:

```text
genomes/genomic/
genomes/cds/
genomes/rna/
```

Reference snapshot: the Borrelia NCBI download used when this guide was
validated contained 82 FASTA files in each folder. The archive is generated on
demand, is not committed to the repository, and may change as NCBI data
changes. The counts below are therefore snapshot expectations, not permanent
invariants.

The CI workflow caches this archive by OS and configuration. On a cache miss,
GitHub Actions downloads it from NCBI before running the test:

```bash
pixi run download-ci-test-data
pixi run ci
```

`functional-test-ci` only consumes an archive that has already been prepared;
it does not download data itself:

```bash
pixi run functional-test-ci
```

## Step 2. Run Through The API

This prepares the dataset, runs the workflow, and verifies the outputs:

```bash
pixi run functional-test
```

Equivalent direct CLI:

```bash
pixi run python -m amprime functional-test
```

After installing the conda package, use:

```bash
amprime functional-test
```

The functional test checks that the workflow and reports complete. Candidate
counts are reported for inspection but are not fixed, because CI may download
a newer NCBI dataset:

```text
functional test ok
primer_rows=<non-negative integer> pcr_rows=<non-negative integer> backend=python report_bytes=<positive integer>
```

## Step 3. Use The Python API

```python
from amprime import AmPrimeProject

project = AmPrimeProject()
result = project.run_functional_test()
print(result.report_html)
```

## Step 4. Manual Commands

Use these only when you want to debug each phase separately.

Run from the repository root.

PowerShell:

```powershell
New-Item -ItemType Directory -Force -Path results\Borrelia | Out-Null
tar -xzf data\borrelia-genomes.tar.gz -C results\Borrelia
```

Bash:

```bash
mkdir -p results/Borrelia
tar -xzf data/borrelia-genomes.tar.gz -C results/Borrelia
```

```bash
pixi run snakemake --cores 4 --rerun-incomplete \
  results/Borrelia/reports/recG_report.html \
  results/Borrelia/reports/gene_report_cross.html
```

Expected rules:

```text
gene_extract -> cluster -> align -> primers_design -> primers_check -> in_silico_pcr -> gene_report -> gene_report_cross
```

`genomes_download` should not run.

## Step 5. Verify Outputs Manually

```bash
pixi run python -m amprime verify --genus Borrelia --gene recG
```

For the reference snapshot only, add `--expect-no-candidates` to turn the
no-candidate result into an explicit snapshot assertion:

```bash
pixi run python -m amprime verify --genus Borrelia --gene recG --expect-no-candidates
```

## Step 6. Interpret The Result

For the reference snapshot and default config, the diagnostic result is
expected to contain no primer candidates:

- `gene_extract` finds 82 `recG` CDS hits.
- `cluster` reduces them to 12 centroids.
- `align` uses the Python backend.
- `gene_report` writes a complete no-candidate report.

These counts and the no-candidate result are not guaranteed for a fresh NCBI
download. The CI functional test intentionally does not assert them. Synthetic
positive primer/PCR behavior is covered by:

```bash
pixi run smoke
```
