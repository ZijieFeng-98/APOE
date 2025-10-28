#!/usr/bin/env python3
"""Utilities to produce APOE exon slices for PI and scientific validation patches.

This module codifies the step-by-step shell guidance that was previously
provided informally.  It exposes two subcommands:

```
patch-a  # PI demonstration – extract exon 2 reads and reference sequence
patch-b  # Scientific validation – exon 4 extraction + rs429358/rs7412 calling
```

Both subcommands share the same input arguments (BAM/CRAM, reference FASTA,
GTF) and emit reproducible artefacts in an output directory.  External tools
(`samtools`, `bedtools`, `bcftools`) are invoked directly so the behaviour
matches the manual recipe exactly.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

__all__ = ["main"]


@dataclass(frozen=True)
class Exon:
    chrom: str
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    exon_number: str
    transcript_id: str

    @property
    def bed_interval(self) -> Tuple[str, int, int]:
        """Return a BED-style (chrom, start-1, end) tuple."""

        return (self.chrom, self.start - 1, self.end)

    @property
    def region_string(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"


@dataclass(frozen=True)
class SNPPosition:
    chrom: str
    position: int  # 1-based genomic coordinate
    rsid: str

    @property
    def bed_row(self) -> Tuple[str, int, int, str]:
        return (self.chrom, self.position - 1, self.position, self.rsid)


BUILD_TO_SNPS: Dict[str, Tuple[SNPPosition, SNPPosition]] = {
    "hg38": (
        SNPPosition("chr19", 44908684, "rs429358"),
        SNPPosition("chr19", 44908822, "rs7412"),
    ),
    "hg19": (
        SNPPosition("chr19", 45411941, "rs429358"),
        SNPPosition("chr19", 45412079, "rs7412"),
    ),
}


class PatchError(RuntimeError):
    """Raised when a required resource is missing or invalid."""


class CommandRunner:
    """Wrapper around subprocess.run with logging, dry-run, and error handling."""

    def __init__(self, dry_run: bool = False) -> None:
        self.dry_run = dry_run

    def run(
        self,
        command: Sequence[str],
        *,
        capture_output: bool = False,
        text: bool = False,
        stdout=None,
        stdin=None,
        input: Optional[str] = None,
    ) -> subprocess.CompletedProcess:
        if self.dry_run:
            print("DRY-RUN:", " ".join(command))
            return subprocess.CompletedProcess(command, 0, "" if text else b"", "" if text else b"")
        try:
            return subprocess.run(
                command,
                check=True,
                capture_output=capture_output,
                text=text,
                stdout=stdout,
                stdin=stdin,
                input=input,
            )
        except subprocess.CalledProcessError as exc:  # pragma: no cover - delegated command
            stderr = exc.stderr if text else (exc.stderr.decode() if exc.stderr else "")
            raise PatchError(f"Command failed ({exc.returncode}): {' '.join(command)}\n{stderr}") from exc

    def pipeline(self, commands: Sequence[Sequence[str]], output_path: Optional[Path] = None) -> None:
        if self.dry_run:
            print("DRY-RUN pipeline:")
            for command in commands:
                print("  ", " ".join(command))
            return

        processes = []
        prev_stdout = None
        try:
            for i, command in enumerate(commands):
                stdout = subprocess.PIPE if i < len(commands) - 1 or output_path else None
                if i == len(commands) - 1 and output_path:
                    output_handle = output_path.open("wb")
                    stdout = output_handle
                else:
                    output_handle = None
                proc = subprocess.Popen(command, stdin=prev_stdout, stdout=stdout, stderr=subprocess.PIPE)
                if prev_stdout is not None and prev_stdout is not subprocess.PIPE:
                    prev_stdout.close()
                prev_stdout = proc.stdout
                processes.append((proc, output_handle))

            for proc, output_handle in processes:
                _, stderr_data = proc.communicate()
                if output_handle:
                    output_handle.close()
                if proc.returncode != 0:
                    raise PatchError(
                        f"Command failed ({proc.returncode}): {' '.join(proc.args)}\n{(stderr_data or b'').decode()}"
                    )
        finally:
            for proc, output_handle in processes:
                if output_handle and not output_handle.closed:
                    output_handle.close()

    __call__ = run


def ensure_tools_available(tools: Iterable[str]) -> None:
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        raise PatchError(f"Missing required tools: {', '.join(missing)}")


def load_apoe_exons(gtf_path: Path, transcript_id: Optional[str] = None) -> Dict[str, Exon]:
    """Parse GTF file and return exon definitions for APOE."""

    transcripts: Dict[str, Dict[str, Exon]] = {}
    with gtf_path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            chrom, _source, feature, start_str, end_str, _score, _strand, _frame, attributes = row
            if feature != "exon":
                continue
            if "APOE" not in attributes:
                continue

            attrs = parse_gtf_attributes(attributes)
            exon_number = attrs.get("exon_number")
            transcript = attrs.get("transcript_id")
            if exon_number is None or transcript is None:
                continue
            if transcript_id and transcript != transcript_id:
                continue

            start = int(start_str)
            end = int(end_str)
            transcripts.setdefault(transcript, {})[exon_number] = Exon(
                chrom=chrom,
                start=start,
                end=end,
                exon_number=exon_number,
                transcript_id=transcript,
            )

    if not transcripts:
        raise PatchError(
            "Unable to locate APOE exons in the supplied GTF. Provide a file aligned with the reference build and "
            "optionally specify --transcript-id if the canonical transcript is named differently."
        )

    if transcript_id:
        try:
            return transcripts[transcript_id]
        except KeyError as exc:
            raise PatchError(
                f"Transcript '{transcript_id}' not found in {gtf_path}. Check the identifier or omit --transcript-id."
            ) from exc

    # Prefer transcripts that contain exon 4 (rs429358/rs7412) and have the most exons annotated.
    canonical_transcript = None
    for candidate, exons in transcripts.items():
        if canonical_transcript is None:
            canonical_transcript = candidate
            continue
        current = transcripts[canonical_transcript]
        candidate_score = ("4" in exons, len(exons))
        current_score = ("4" in current, len(current))
        if candidate_score > current_score:
            canonical_transcript = candidate

    return transcripts[canonical_transcript]


def parse_gtf_attributes(attributes: str) -> Dict[str, str]:
    entries = {}
    for field in attributes.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        if " " not in field:
            continue
        key, value = field.split(" ", 1)
        entries[key] = value.strip().strip('"')
    return entries


def write_bed(exon: Exon, path: Path, name: str) -> Path:
    content = "\t".join((*map(str, exon.bed_interval), name)) + "\n"
    path.write_text(content)
    return path


def write_snp_bed(snps: Sequence[SNPPosition], path: Path) -> Path:
    lines = ["\t".join(str(value) for value in snp.bed_row) for snp in snps]
    path.write_text("\n".join(lines) + "\n")
    return path


def slice_bam_to_file(runner: CommandRunner, bam_path: Path, bed_path: Path, output_bam: Path) -> None:
    with output_bam.open("wb") as handle:
        runner.run(
            ["samtools", "view", "-b", "-L", str(bed_path), str(bam_path)],
            stdout=handle,
        )
    runner.run(["samtools", "index", str(output_bam)])


def extract_reference_sequence(runner: CommandRunner, ref_fasta: Path, bed_path: Path, output_fasta: Path) -> None:
    runner.run(
        ["bedtools", "getfasta", "-fi", str(ref_fasta), "-bed", str(bed_path), "-name", "-fo", str(output_fasta)]
    )


def compute_depth(runner: CommandRunner, bam_path: Path, region: str) -> float:
    proc = runner.run(["samtools", "depth", "-aa", "-r", region, str(bam_path)], capture_output=True, text=True)
    if not proc.stdout:
        return 0.0
    total = 0
    count = 0
    for line in proc.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 3:
            total += int(parts[2])
            count += 1
    return total / count if count else 0.0


def call_two_snps(
    runner: CommandRunner,
    bam_path: Path,
    ref_fasta: Path,
    snp_bed: Path,
    vcf_path: Path,
    tsv_path: Path,
) -> None:
    pipeline_commands = [
        [
            "bcftools",
            "mpileup",
            "-f",
            str(ref_fasta),
            "-R",
            str(snp_bed),
            str(bam_path),
        ],
        ["bcftools", "call", "-mv"],
        ["bcftools", "view", "-Oz"],
    ]
    runner.pipeline(pipeline_commands, output_path=vcf_path)
    runner.run(["bcftools", "index", str(vcf_path)])
    query = runner.run(
        [
            "bcftools",
            "query",
            "-f",
            "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%DP\t%GQ\t%GT]\n",
            str(vcf_path),
        ],
        capture_output=True,
        text=True,
    )
    tsv_path.write_text(query.stdout)


def run_patch(
    patch: str,
    bam_paths: Sequence[Path],
    labels: Sequence[str],
    ref_fasta: Path,
    gtf_path: Path,
    output_dir: Path,
    transcript_id: Optional[str],
    build: str,
    dry_run: bool = False,
) -> Dict[str, Dict[str, str]]:
    ensure_tools_available(["samtools", "bedtools"] + (["bcftools"] if patch == "patch-b" else []))

    exons = load_apoe_exons(gtf_path, transcript_id=transcript_id)
    try:
        exon2 = exons["2"]
    except KeyError as exc:
        raise PatchError("The APOE exon 2 entry was not found in the provided GTF.") from exc
    try:
        exon4 = exons["4"]
    except KeyError:
        exon4 = None

    runner = CommandRunner(dry_run=dry_run)

    output_dir.mkdir(parents=True, exist_ok=True)

    results: Dict[str, Dict[str, str]] = {}

    for bam_path, label in zip(bam_paths, labels):
        sample_dir = output_dir / label
        sample_dir.mkdir(exist_ok=True)
        artifacts: Dict[str, str] = {}

        # Patch A – always perform exon 2 extraction
        exon2_bed = sample_dir / f"{label}_exon2.bed"
        write_bed(exon2, exon2_bed, "APOE_exon2")

        exon2_bam = sample_dir / f"{label}_APOE_exon2.bam"
        slice_bam_to_file(runner, bam_path, exon2_bed, exon2_bam)
        artifacts["exon2_bam"] = str(exon2_bam)

        exon2_fasta = sample_dir / f"{label}_APOE_exon2.fasta"
        extract_reference_sequence(runner, ref_fasta, exon2_bed, exon2_fasta)
        artifacts["exon2_fasta"] = str(exon2_fasta)

        depth = compute_depth(runner, exon2_bam, exon2.region_string)
        artifacts["exon2_mean_depth"] = f"{depth:.2f}"

        if patch == "patch-a":
            results[label] = artifacts
            continue

        if patch == "patch-b":
            if exon4 is None:
                raise PatchError("The GTF does not contain exon 4 for APOE; cannot run patch-b.")

            exon4_bed = sample_dir / f"{label}_exon4.bed"
            write_bed(exon4, exon4_bed, "APOE_exon4")

            exon4_bam = sample_dir / f"{label}_APOE_exon4.bam"
            slice_bam_to_file(runner, bam_path, exon4_bed, exon4_bam)
            artifacts["exon4_bam"] = str(exon4_bam)

            exon4_fasta = sample_dir / f"{label}_APOE_exon4.fasta"
            extract_reference_sequence(runner, ref_fasta, exon4_bed, exon4_fasta)
            artifacts["exon4_fasta"] = str(exon4_fasta)

            exon4_depth = compute_depth(runner, exon4_bam, exon4.region_string)
            artifacts["exon4_mean_depth"] = f"{exon4_depth:.2f}"

            snps = BUILD_TO_SNPS.get(build)
            if not snps:
                raise PatchError(f"Unknown build '{build}'. Expected one of: {', '.join(BUILD_TO_SNPS)}")
            snp_bed = sample_dir / f"{label}_apoe_two_snps.bed"
            write_snp_bed(snps, snp_bed)

            vcf_path = sample_dir / f"{label}_apoe_two_snps.vcf.gz"
            tsv_path = sample_dir / f"{label}_apoe_two_snps.tsv"
            call_two_snps(runner, exon4_bam, ref_fasta, snp_bed, vcf_path, tsv_path)
            artifacts["apoe_two_snps_vcf"] = str(vcf_path)
            artifacts["apoe_two_snps_tsv"] = str(tsv_path)

            results[label] = artifacts
        else:
            raise PatchError(f"Unsupported patch '{patch}'.")

    summary_path = output_dir / f"{patch}_summary.json"
    summary_path.write_text(json.dumps(results, indent=2))
    return results


def normalise_labels(bam_paths: Sequence[Path], labels: Optional[Sequence[str]]) -> List[str]:
    if labels:
        if len(labels) != len(bam_paths):
            raise PatchError("The number of --labels must match the number of --bam entries.")
        return list(labels)
    return [bam_path.stem for bam_path in bam_paths]


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="patch", required=True)

    def add_common_arguments(subparser: argparse.ArgumentParser) -> None:
        subparser.add_argument("--bam", action="append", required=True, type=Path, help="Input BAM/CRAM file(s).")
        subparser.add_argument(
            "--label",
            action="append",
            dest="labels",
            help="Sample labels corresponding to the --bam files (defaults to BAM stem).",
        )
        subparser.add_argument("--reference", required=True, type=Path, help="Reference FASTA file.")
        subparser.add_argument("--gtf", required=True, type=Path, help="Annotation GTF file containing APOE exons.")
        subparser.add_argument(
            "--output-dir",
            default=Path("apoe_patch_outputs"),
            type=Path,
            help="Directory for generated artefacts.",
        )
        subparser.add_argument(
            "--transcript-id",
            help="Optional transcript identifier to pin APOE exon selection (e.g., ENST00000252486).",
        )
        subparser.add_argument(
            "--build",
            default="hg38",
            choices=sorted(BUILD_TO_SNPS.keys()),
            help="Reference build to determine SNP coordinates (patch-b only).",
        )
        subparser.add_argument("--dry-run", action="store_true", help="Print commands without executing them.")

    add_common_arguments(subparsers.add_parser("patch-a", help="PI demonstration (exon 2 extraction)."))
    add_common_arguments(subparsers.add_parser("patch-b", help="Scientific validation (exon 4 + SNP calling)."))
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)

    bam_paths = [path.resolve() for path in args.bam]
    labels = normalise_labels(bam_paths, args.labels)

    try:
        run_patch(
            patch=args.patch,
            bam_paths=bam_paths,
            labels=labels,
            ref_fasta=args.reference.resolve(),
            gtf_path=args.gtf.resolve(),
            output_dir=args.output_dir.resolve(),
            transcript_id=args.transcript_id,
            build=args.build,
            dry_run=args.dry_run,
        )
    except PatchError as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

