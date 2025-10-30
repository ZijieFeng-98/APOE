#!/usr/bin/env python3
"""Generate a clinical-style APOE interpretation from VCF calls.

The original version of this script assumed that two separate VCF files were
provided—one for rs429358 and one for rs7412—and attempted to derive the APOE
genotype by pairing alleles by index.  That approach has a couple of major
problems:

* It silently accepts bcftools error files (which look like VCFs but contain no
  data) and reports "Insufficient data" without explaining the failure.
* It pairs the first allele of rs429358 with the first allele of rs7412 even
  though unphased genotypes provide no such linkage.  This led to incorrect
  allele assignments (e.g., ε4 being reported as ε2).

This rewritten module provides a small but more resilient APOE interpretation
engine that understands the biology of the three alleles (ε2, ε3, ε4) and
counts them directly from the underlying SNP genotypes.  The CLI now accepts a
single VCF or the two-file workflow, emits detailed warnings in the report, and
remains backwards compatible with the legacy positional arguments used in the
existing pipeline scripts.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

__all__ = ["main"]


class InterpretationError(RuntimeError):
    """Raised when the input data cannot support an APOE interpretation."""


@dataclass(frozen=True)
class VariantTarget:
    """Define a variant of interest."""

    name: str
    rsid: str
    position: int


@dataclass
class VariantCall:
    """Minimal representation of a single-sample VCF call."""

    target: VariantTarget
    chromosome: str
    position: int
    rsid: Optional[str]
    ref: str
    alts: List[str]
    genotype: str
    alleles: List[str]
    sample_id: Optional[str]


APOE_TARGETS: Tuple[VariantTarget, VariantTarget] = (
    VariantTarget(name="rs429358", rsid="rs429358", position=45411941),
    VariantTarget(name="rs7412", rsid="rs7412", position=45412079),
)


def _extract_sample_id(header_line: str) -> Optional[str]:
    columns = header_line.rstrip().split("\t")
    if len(columns) > 9:
        return columns[9]
    return None


def _decode_genotype(ref: str, alt_field: str, genotype: str) -> Optional[List[str]]:
    if not genotype or genotype == "./.":
        return None

    alleles: List[str] = []
    alt_alleles = alt_field.split(",") if alt_field else []
    for code in re.split(r"[\/|]", genotype):
        if code == ".":
            return None
        try:
            index = int(code)
        except ValueError:
            return None
        if index == 0:
            alleles.append(ref)
        else:
            alt_index = index - 1
            if alt_index >= len(alt_alleles):
                return None
            alleles.append(alt_alleles[alt_index])
    return alleles


def _target_lookup(targets: Iterable[VariantTarget]) -> Tuple[Dict[str, VariantTarget], Dict[int, VariantTarget]]:
    by_rsid: Dict[str, VariantTarget] = {}
    by_position: Dict[int, VariantTarget] = {}
    for target in targets:
        by_rsid[target.rsid] = target
        by_position[target.position] = target
    return by_rsid, by_position


def read_targets_from_vcf(vcf_path: Path, targets: Iterable[VariantTarget]) -> Dict[str, VariantCall]:
    """Parse a VCF and return calls for the requested targets."""

    calls: Dict[str, VariantCall] = {}
    by_rsid, by_position = _target_lookup(targets)
    sample_id: Optional[str] = None

    try:
        content = vcf_path.read_text().splitlines()
    except OSError as exc:  # pragma: no cover - defensive IO path
        raise InterpretationError(f"Unable to read VCF '{vcf_path}': {exc}") from exc

    if any(line.startswith("[E::hts_open_format]") for line in content):
        raise InterpretationError(
            f"File '{vcf_path}' appears to contain an error message from bcftools instead of VCF data."
        )

    for raw_line in content:
        if raw_line.startswith("##"):
            continue
        if raw_line.startswith("#CHROM"):
            sample_id = _extract_sample_id(raw_line)
            continue
        if not raw_line.strip():
                    continue

        fields = raw_line.split("\t")
                if len(fields) < 10:
            # Not a complete single-sample VCF record.
            continue

        chrom, pos_str, rsid, ref, alt, _qual, _filter, _info, fmt, sample = fields[:10]
        try:
            position = int(pos_str)
        except ValueError:
            continue

        target = by_rsid.get(rsid)
        if not target:
            target = by_position.get(position)
        if not target:
                    continue

        format_fields = fmt.split(":")
        sample_fields = sample.split(":")
        try:
            gt_index = format_fields.index("GT")
        except ValueError:
            raise InterpretationError(
                f"Record for {target.name} in '{vcf_path}' is missing a GT field."
            )

        try:
                    genotype = sample_fields[gt_index]
        except IndexError:
            raise InterpretationError(
                f"Sample column for {target.name} in '{vcf_path}' does not contain GT values."
            )

        alleles = _decode_genotype(ref, alt, genotype)
        if not alleles:
            raise InterpretationError(
                f"Could not decode alleles for {target.name} in '{vcf_path}' (genotype '{genotype}')."
            )

        calls[target.name] = VariantCall(
            target=target,
            chromosome=chrom,
            position=position,
            rsid=None if rsid == "." else rsid,
            ref=ref,
            alts=alt.split(",") if alt else [],
            genotype=genotype,
            alleles=alleles,
            sample_id=sample_id,
        )

    return calls


def determine_apoe(call_429358: VariantCall, call_7412: VariantCall) -> Tuple[str, List[str]]:
    """Determine the APOE diplotype and return it with interpretation notes."""

    alleles_429358 = call_429358.alleles
    alleles_7412 = call_7412.alleles

    if len(alleles_429358) != 2 or len(alleles_7412) != 2:
        raise InterpretationError("Expected diploid genotypes for both variants.")

    e2_count = alleles_7412.count("T")
    e4_count = alleles_429358.count("C")
    e3_count = 2 - e2_count - e4_count

    if e3_count < 0 or e3_count > 2:
        raise InterpretationError(
            "Observed allele counts are incompatible with APOE haplotypes. Check for sample mix-up or low-quality calls."
        )

    genotype_parts = ["ε2"] * e2_count + ["ε4"] * e4_count + ["ε3"] * e3_count
    genotype_parts.sort()

    if len(genotype_parts) != 2:
        raise InterpretationError("Unable to resolve a diploid APOE genotype from the provided data.")

    notes: List[str] = []
    if e2_count and e4_count:
        notes.append(
            "Both ε2 and ε4 alleles detected; consider confirming phasing if clinically indicated."
        )

    return "/".join(genotype_parts), notes


def clinical_interp(genotype: str) -> Tuple[str, str, str]:
    interp = {
        "ε2/ε2": ("REDUCED RISK", "Protective relative to population baseline", "~0.5x"),
        "ε2/ε3": ("REDUCED RISK", "Somewhat protective", "~0.6x"),
        "ε2/ε4": ("MODERATE RISK", "Mixed effect; phenotype can depend on other factors", "~2-3x"),
        "ε3/ε3": ("AVERAGE RISK", "Most common APOE genotype", "1x"),
        "ε3/ε4": ("INCREASED RISK", "Moderately increased risk", "~3x"),
        "ε4/ε4": ("HIGH RISK", "Significantly increased population risk", "~8-12x"),
    }
    return interp.get(genotype, ("UNKNOWN", "Unable to determine", "Unknown"))


def _format_variant_block(call: VariantCall) -> str:
    alleles = "/".join(call.alleles)
    rsid = call.rsid or call.target.rsid
    return (
        f"{call.target.name} ({rsid}; chr{call.chromosome}:{call.position}):\n"
        f"  Reference Allele: {call.ref}\n"
        f"  Alternate Alleles: {', '.join(call.alts) if call.alts else '—'}\n"
        f"  Genotype: {alleles} ({call.genotype})\n"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Interpret APOE genotype from small VCF extracts produced by whole-genome sequencing."
        )
    )
    parser.add_argument(
        "--vcf",
        type=Path,
        help=(
            "Single VCF file containing both rs429358 and rs7412 (preferred). "
            "If omitted, provide individual files with --rs429358 and --rs7412."
        ),
    )
    parser.add_argument("--rs429358", type=Path, help="VCF containing the rs429358 call")
    parser.add_argument("--rs7412", type=Path, help="VCF containing the rs7412 call")
    parser.add_argument("-o", "--output", type=Path, help="Path to write the report")
    parser.add_argument("--sample-id", help="Sample identifier to include in the report")
    parser.add_argument(
        "--analysis-date",
        help="Override the analysis date in YYYY-MM-DD format (defaults to today)",
    )
    parser.add_argument(
        "legacy_args",
        nargs="*",
        help=argparse.SUPPRESS,
    )
    return parser


def _coerce_date(value: Optional[str]) -> str:
    if not value:
        return _dt.date.today().isoformat()
    try:
        return _dt.date.fromisoformat(value).isoformat()
    except ValueError as exc:
        raise InterpretationError(
            f"Analysis date '{value}' is not a valid ISO date (YYYY-MM-DD)."
        ) from exc


def _resolve_args(args: argparse.Namespace) -> Tuple[Path, Path, Optional[Path], Optional[str], str]:
    # Legacy positional mode: interpret.py rs429358.vcf rs7412.vcf output [date]
    if args.legacy_args and not any((args.vcf, args.rs429358, args.rs7412, args.output)):
        if len(args.legacy_args) not in {3, 4}:
            raise InterpretationError(
                "Legacy invocation expects 3 or 4 positional arguments: rs429358.vcf rs7412.vcf output [date]"
            )
        rs429358_path = Path(args.legacy_args[0])
        rs7412_path = Path(args.legacy_args[1])
        output_path = Path(args.legacy_args[2])
        analysis_date = _coerce_date(args.legacy_args[3] if len(args.legacy_args) == 4 else None)
        return rs429358_path, rs7412_path, output_path, None, analysis_date

    output_path = args.output
    if output_path is None:
        raise InterpretationError("An output path is required (use --output or legacy positional argument).")

    analysis_date = _coerce_date(args.analysis_date)

    if args.vcf:
        return args.vcf, args.vcf, output_path, args.sample_id, analysis_date

    if not (args.rs429358 and args.rs7412):
        raise InterpretationError("Provide both --rs429358 and --rs7412 when --vcf is not supplied.")

    return args.rs429358, args.rs7412, output_path, args.sample_id, analysis_date


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    try:
        rs429358_path, rs7412_path, output_path, sample_override, analysis_date = _resolve_args(args)
    except InterpretationError as exc:
        parser.error(str(exc))

    collected_calls: Dict[str, VariantCall] = {}
    errors: List[str] = []
    sample_id: Optional[str] = None

    for path, target in (
        (rs429358_path, (APOE_TARGETS[0],)),
        (rs7412_path, (APOE_TARGETS[1],)),
    ):
        try:
            calls = read_targets_from_vcf(Path(path), target)
        except InterpretationError as exc:
            errors.append(str(exc))
            continue
        collected_calls.update(calls)
        sample_id = sample_id or next(iter(calls.values())).sample_id

    if APOE_TARGETS[0].name not in collected_calls or APOE_TARGETS[1].name not in collected_calls:
        missing = [
            t.name
            for t in APOE_TARGETS
            if t.name not in collected_calls
        ]
        errors.append(f"Missing required variant calls: {', '.join(missing)}")

    if errors:
        # When using a single VCF, both targets come from the same file and we only
        # want to report the failure once.
        error_text = "\n".join(errors)
        raise SystemExit(f"Unable to interpret APOE genotype:\n{error_text}")

    rs429358_call = collected_calls[APOE_TARGETS[0].name]
    rs7412_call = collected_calls[APOE_TARGETS[1].name]

    if sample_override:
        sample_id = sample_override

    try:
        genotype, genotype_notes = determine_apoe(rs429358_call, rs7412_call)
    except InterpretationError as exc:
        raise SystemExit(f"Unable to interpret APOE genotype:\n{exc}") from exc

risk, desc, rel_risk = clinical_interp(genotype)

    output_lines = ["=" * 80, "APOE GENOTYPING REPORT", "=" * 80, ""]
    output_lines.append(f"Sample ID: {sample_id or 'Unknown Sample'}")
    output_lines.append(f"Analysis Date: {analysis_date}")
    output_lines.append("Reference: GRCh37/hg19")
    output_lines.extend(["", "-" * 80, "RAW GENOTYPE DATA", "-" * 80, ""])
    output_lines.append(_format_variant_block(rs429358_call))
    output_lines.append(_format_variant_block(rs7412_call))
    output_lines.extend(["-" * 80, "APOE GENOTYPE", "-" * 80, "", f"APOE Genotype: {genotype}", ""])
    output_lines.extend(["-" * 80, "CLINICAL INTERPRETATION", "-" * 80, ""])
    output_lines.append(f"Risk Category: {risk}")
    output_lines.append(f"Description: {desc}")
    output_lines.append(f"Relative Risk: {rel_risk}")
    output_lines.append("")

    if genotype_notes:
        output_lines.append("Genotype Notes:")
        for note in genotype_notes:
            output_lines.append(f"- {note}")
        output_lines.append("")

    output_lines.append("Important Notes:")
    output_lines.append("- APOE is just one of many genetic and lifestyle risk factors")
    output_lines.append("- This interpretation is for research purposes only")
    output_lines.append("- Consult a clinical geneticist or genetic counselor for medical decisions")
    output_lines.append("")
    output_lines.append("=" * 80)

    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text("\n".join(output_lines) + "\n")

print(f"APOE Genotype: {genotype} | Risk: {risk}")


if __name__ == "__main__":
    main()
