# Root entry point for Snakemake.
# The workflow implementation lives under workflow/ so project files stay tidy,
# but users can still run `snakemake -n` from the repository root.

include: "workflow/Snakefile"
