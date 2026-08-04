"""Public Python API for AmPrime workflow orchestration and result checks.

The API intentionally stays thin: it prepares local input data, calls
Snakemake as the scheduler, and inspects the files produced by the existing
workflow scripts.
"""

from __future__ import annotations

import contextlib
import csv
import os
import shutil
import subprocess
import tarfile
from dataclasses import dataclass
from datetime import UTC, datetime
from importlib import resources
from pathlib import Path

from .provenance import fasta_directory_summary, sha256_file

DEFAULT_GENUS = "Borrelia"
DEFAULT_GENE = "recG"
DEFAULT_TEST_ARCHIVE = Path("data") / "borrelia-genomes.tar.gz"


@dataclass(frozen=True)
class ResultPaths:
    """Common per-gene result paths."""

    report_html: Path
    primers_tsv: Path
    amplicons_tsv: Path
    species_summary_tsv: Path
    species_tsv: Path
    diversity_png: Path
    alignment_tsv: Path


@dataclass(frozen=True)
class PipelineRun:
    """Completed Snakemake invocation."""

    command: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str
    target: str | None


@dataclass(frozen=True)
class FunctionalTestResult:
    """Summary of the local functional test output."""

    genus: str
    gene: str
    report_html: Path
    report_bytes: int
    primer_rows: int
    pcr_rows: int
    alignment_backend: str
    requested_backend: str


def _project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _resource_path(*parts: str) -> Path:
    source_path = _project_root().joinpath(*parts)
    if source_path.exists():
        return source_path
    return Path(str(resources.files("amprime").joinpath(*parts)))


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _safe_extract_tar_gz(archive: Path, destination: Path) -> None:
    destination = destination.resolve()
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar.getmembers():
            if member.issym() or member.islnk():
                raise ValueError(
                    f"Refusing to extract link from archive: {member.name}"
                )
            if not (member.isdir() or member.isfile()):
                raise ValueError(
                    f"Refusing to extract special archive member: {member.name}"
                )
            target = (destination / member.name).resolve()
            if target != destination and destination not in target.parents:
                raise ValueError(
                    f"Refusing to extract unsafe archive member: {member.name}"
                )
        tar.extractall(destination, filter="data")


def _write_local_manifest(genomes_dir: Path, genus: str, config: Path) -> None:
    generated_at = datetime.now(UTC).isoformat(timespec="seconds")
    rows = []
    for label, fmt in [
        ("genomic", "fasta"),
        ("cds", "cds-fasta"),
        ("rna", "rna-fasta"),
    ]:
        directory = genomes_dir / label
        rows.append(
            {
                "label": label,
                "format": fmt,
                "output_dir": str(directory),
                **fasta_directory_summary(directory),
            }
        )

    manifest = genomes_dir / "download_manifest.tsv"
    fieldnames = [
        "generated_at",
        "genus",
        "assembly_level",
        "label",
        "format",
        "output_dir",
        "n_fna",
        "total_bytes",
        "data_fingerprint",
        "config_sha256",
    ]
    with manifest.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "generated_at": generated_at,
                    "genus": genus,
                    "assembly_level": "local-archive",
                    "config_sha256": sha256_file(config),
                    **row,
                }
            )


def _as_snakemake_target(path: str | Path | None, root: Path) -> str | None:
    if path is None:
        return None
    target = Path(path)
    if target.is_absolute():
        with contextlib.suppress(ValueError):
            target = target.relative_to(root)
    return target.as_posix()


class AmPrimeProject:
    """Convenience wrapper around an AmPrime checkout."""

    def __init__(self, root: str | Path | None = None):
        self.root = Path(root).resolve() if root else Path.cwd().resolve()

    def snakefile(self) -> Path:
        root_snakefile = self.root / "Snakefile"
        if root_snakefile.is_file():
            return root_snakefile
        return _resource_path("workflow", "Snakefile")

    def workflow_scripts_dir(self) -> Path:
        candidates = [
            self.root / "workflow" / "scripts",
            self.snakefile().parent / "scripts",
            _resource_path("workflow", "scripts"),
        ]
        for candidate in candidates:
            if (candidate / "config_schema.py").is_file():
                return candidate
        searched = "\n".join(str(candidate) for candidate in candidates)
        raise ModuleNotFoundError(
            f"Cannot find config_schema.py. Searched:\n{searched}"
        )

    def ensure_default_config(self) -> Path:
        config_path = self.root / "config" / "config.yaml"
        if config_path.is_file():
            return config_path
        source_config = _resource_path("config", "config.yaml")
        config_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_config, config_path)
        return config_path

    def result_dir(self, genus: str = DEFAULT_GENUS) -> Path:
        return self.root / "results" / genus

    def result_paths(
        self, genus: str = DEFAULT_GENUS, gene: str = DEFAULT_GENE
    ) -> ResultPaths:
        result_dir = self.result_dir(genus)
        return ResultPaths(
            report_html=result_dir / "reports" / f"{gene}_report.html",
            primers_tsv=result_dir / "primers" / f"{gene}_primers.tsv",
            amplicons_tsv=result_dir / "primers" / f"{gene}_amplicons.tsv",
            species_summary_tsv=result_dir / "primers" / f"{gene}_species_summary.tsv",
            species_tsv=result_dir / "primers" / f"{gene}_species.tsv",
            diversity_png=result_dir / "primers" / f"{gene}_diversity.png",
            alignment_tsv=result_dir / "aligned" / f"{gene}.alignment.tsv",
        )

    def prepare_local_dataset(
        self, archive: str | Path = DEFAULT_TEST_ARCHIVE, genus: str = DEFAULT_GENUS
    ) -> Path:
        archive_path = Path(archive)
        if not archive_path.is_absolute():
            archive_path = self.root / archive_path
        if not archive_path.is_file():
            raise FileNotFoundError(f"Cannot find test dataset archive: {archive_path}")

        config_path = self.ensure_default_config()
        output_dir = self.result_dir(genus)
        output_dir.mkdir(parents=True, exist_ok=True)
        genomes_dir = output_dir / "genomes"
        if genomes_dir.exists():
            shutil.rmtree(genomes_dir)
        _safe_extract_tar_gz(archive_path, output_dir)

        for name in ["genomic", "cds", "rna"]:
            folder = genomes_dir / name
            if not folder.is_dir():
                raise FileNotFoundError(f"Dataset is missing {folder}")
            if not list(folder.glob("*.fna")):
                raise FileNotFoundError(f"Dataset has no FASTA files in {folder}")
        _write_local_manifest(genomes_dir, genus, config_path)
        return genomes_dir

    def run_pipeline(
        self,
        target: str | Path | None = None,
        cores: int = 4,
        dry_run: bool = False,
        rerun_incomplete: bool = True,
        extra_args: list[str] | None = None,
    ) -> PipelineRun:
        self.ensure_default_config()
        command = ["snakemake", "-s", str(self.snakefile())]
        if dry_run:
            command.append("-n")
        else:
            command.extend(["--cores", str(cores)])
            if rerun_incomplete:
                command.append("--rerun-incomplete")
        if extra_args:
            command.extend(extra_args)
        snakemake_target = _as_snakemake_target(target, self.root)
        if snakemake_target is not None:
            command.append(snakemake_target)

        env = os.environ.copy()
        scripts_dir = str(self.workflow_scripts_dir())
        pythonpath = env.get("PYTHONPATH")
        env["PYTHONPATH"] = (
            scripts_dir if not pythonpath else f"{scripts_dir}{os.pathsep}{pythonpath}"
        )

        completed = subprocess.run(  # noqa: S603 - executable is fixed to Snakemake; shell is disabled.
            command, cwd=self.root, check=False, text=True, capture_output=True, env=env
        )
        if completed.returncode != 0:
            raise subprocess.CalledProcessError(
                completed.returncode,
                command,
                output=completed.stdout,
                stderr=completed.stderr,
            )
        return PipelineRun(
            command=tuple(command),
            returncode=completed.returncode,
            stdout=completed.stdout,
            stderr=completed.stderr,
            target=snakemake_target,
        )

    def verify_result_outputs(
        self,
        genus: str = DEFAULT_GENUS,
        gene: str = DEFAULT_GENE,
        expect_no_candidates: bool = False,
    ) -> FunctionalTestResult:
        paths = self.result_paths(genus, gene)
        required = [
            paths.report_html,
            paths.primers_tsv,
            paths.amplicons_tsv,
            paths.species_summary_tsv,
            paths.species_tsv,
            paths.diversity_png,
            paths.alignment_tsv,
        ]
        for path in required:
            if not path.is_file():
                raise FileNotFoundError(f"Missing expected output: {path}")
        if paths.diversity_png.stat().st_size <= 0:
            raise AssertionError(f"Empty diversity plot: {paths.diversity_png}")

        cross_gene_report = paths.report_html.parent / "cross_gene_report.html"
        if not cross_gene_report.is_file():
            raise FileNotFoundError(f"Missing expected output: {cross_gene_report}")

        primer_rows = _read_tsv(paths.primers_tsv)
        pcr_rows = _read_tsv(paths.amplicons_tsv)
        alignment_rows = _read_tsv(paths.alignment_tsv)
        if not alignment_rows:
            raise AssertionError(f"Empty alignment metadata: {paths.alignment_tsv}")

        html = paths.report_html.read_text(encoding="utf-8")
        if "Alignment" not in html:
            raise AssertionError("Report is missing the alignment metadata section")
        if "Provenance" not in html or "Config SHA-256" not in html:
            raise AssertionError("Report is missing provenance fingerprints")
        if "Species-level validation" not in html:
            raise AssertionError("Report is missing species-level validation")
        if expect_no_candidates:
            if primer_rows or pcr_rows:
                raise AssertionError("Expected no primer/PCR rows for this dataset")
            if "No candidate primers were found" not in html:
                raise AssertionError("Report is missing the no-candidate message")

        meta = alignment_rows[0]
        return FunctionalTestResult(
            genus=genus,
            gene=gene,
            report_html=paths.report_html,
            report_bytes=paths.report_html.stat().st_size,
            primer_rows=len(primer_rows),
            pcr_rows=len(pcr_rows),
            alignment_backend=meta.get("backend_used", ""),
            requested_backend=meta.get("requested_backend", ""),
        )

    def run_functional_test(
        self,
        archive: str | Path = DEFAULT_TEST_ARCHIVE,
        genus: str = DEFAULT_GENUS,
        gene: str = DEFAULT_GENE,
        cores: int = 4,
    ) -> FunctionalTestResult:
        self.prepare_local_dataset(archive=archive, genus=genus)
        report_target = self.result_paths(genus, gene).report_html
        self.run_pipeline(
            target=report_target,
            cores=cores,
            rerun_incomplete=True,
            extra_args=[
                "--rerun-triggers",
                "mtime",
                "--forcerun",
                "extract_genes",
                "cluster",
                "align",
                "design_primers",
                "check_primers",
                "in_silico_pcr",
                "gene_report",
            ],
        )
        cross_target = self.result_dir(genus) / "reports" / "cross_gene_report.html"
        self.run_pipeline(
            target=cross_target,
            cores=cores,
            rerun_incomplete=True,
            extra_args=[
                "--rerun-triggers",
                "mtime",
                "--forcerun",
                "cross_gene_report",
            ],
        )
        return self.verify_result_outputs(
            genus=genus, gene=gene, expect_no_candidates=True
        )


def prepare_local_dataset(
    archive: str | Path = DEFAULT_TEST_ARCHIVE,
    genus: str = DEFAULT_GENUS,
    root: str | Path | None = None,
) -> Path:
    return AmPrimeProject(root).prepare_local_dataset(archive=archive, genus=genus)


def run_pipeline(
    target: str | Path | None = None,
    cores: int = 4,
    dry_run: bool = False,
    root: str | Path | None = None,
) -> PipelineRun:
    return AmPrimeProject(root).run_pipeline(
        target=target, cores=cores, dry_run=dry_run
    )


def verify_result_outputs(
    genus: str = DEFAULT_GENUS,
    gene: str = DEFAULT_GENE,
    root: str | Path | None = None,
    expect_no_candidates: bool = False,
) -> FunctionalTestResult:
    return AmPrimeProject(root).verify_result_outputs(
        genus=genus, gene=gene, expect_no_candidates=expect_no_candidates
    )


def run_functional_test(
    archive: str | Path = DEFAULT_TEST_ARCHIVE,
    genus: str = DEFAULT_GENUS,
    gene: str = DEFAULT_GENE,
    cores: int = 4,
    root: str | Path | None = None,
) -> FunctionalTestResult:
    return AmPrimeProject(root).run_functional_test(
        archive=archive, genus=genus, gene=gene, cores=cores
    )
