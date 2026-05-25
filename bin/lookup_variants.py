#!/usr/bin/env python3
"""Look up a list of named variants in a sample VCF and (if absent) the CRAM.

Reads a YAML asset (e.g. assets/curated_risk_variants.yaml or
assets/nutrigenomic_variants.yaml) listing variants by GRCh38 coordinates,
and emits a TSV with one row per variant containing the sample's genotype,
risk-allele count, and the curated interpretation text.

Key correctness feature: when a variant is absent from the sample VCF (which
contains only variant sites for DRAGEN-called Nebula data), the script
falls back to samtools mpileup on the CRAM to distinguish:

  - homozygous reference (depth >= threshold, all reads support ref)
  - no-call due to low coverage

Without that fallback, an absent variant looks identical to a hom-ref call
in a downstream tool, which is wrong for clinically actionable variants
where the difference matters.

Combinations declared in the YAML (e.g. APOE ε2/ε3/ε4 diplotype) are passed
to their handler script after the per-variant lookup completes.

Usage:
    python lookup_variants.py
        --variants assets/curated_risk_variants.yaml
        --sample sample.vcf.gz
        --cram sample.cram
        --reference-fasta Homo_sapiens_assembly38.fasta
        --out 04_curated_risk_variants.tsv
        [--out-combinations 04_curated_combinations.tsv]
        [--min-depth 10]
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path
from typing import Optional

try:
    import pysam
    import yaml
except ImportError as exc:
    print(f"ERROR: missing dependency ({exc}). Install pysam + pyyaml.",
          file=sys.stderr)
    sys.exit(1)


def gt_summary(gt_tuple) -> str:
    if gt_tuple is None or any(a is None for a in gt_tuple):
        return "./."
    return "/".join(str(a) for a in gt_tuple)


def zygosity_class(gt_tuple) -> str:
    if gt_tuple is None or any(a is None for a in gt_tuple):
        return "nocall"
    if len(gt_tuple) == 1:
        return "hemi_alt" if gt_tuple[0] >= 1 else "hemi_ref"
    a, b = gt_tuple[0], gt_tuple[1]
    if a == 0 and b == 0:
        return "hom_ref"
    if a >= 1 and b >= 1:
        return "hom_alt"
    return "het"


def risk_allele_count(zyg: str, gt_tuple, sample_alt_index_is_risk: bool) -> int:
    """Count copies of the risk allele (0, 1, or 2).

    `sample_alt_index_is_risk` is True if the matched alt corresponds to the
    YAML's risk_allele. If risk_allele equals ref, the count is inverted.
    """
    if zyg in ("nocall",):
        return -1
    if zyg == "hom_ref":
        return 0 if sample_alt_index_is_risk else 2
    if zyg == "het":
        return 1
    if zyg == "hom_alt":
        return 2 if sample_alt_index_is_risk else 0
    if zyg == "hemi_alt":
        return 1 if sample_alt_index_is_risk else 0
    if zyg == "hemi_ref":
        return 0 if sample_alt_index_is_risk else 1
    return -1


def cram_pileup_inference(cram: Path, reference: Path, chrom: str, pos: int,
                          ref: str, alt: str, min_depth: int) -> tuple[str, str, int]:
    """When the variant is absent from the VCF, infer call from the CRAM.

    Returns (call_type, raw_pileup, depth). call_type ∈ {hom_ref, low_coverage,
    unexpected_alt, no_data}.
    """
    proc = subprocess.run(
        ["samtools", "mpileup", "-f", str(reference),
         "-r", f"{chrom}:{pos}-{pos}",
         "-q", "20", "-Q", "20", str(cram)],
        capture_output=True, text=True, check=False,
    )
    if proc.returncode != 0 or not proc.stdout.strip():
        return "no_data", "", 0
    parts = proc.stdout.strip().split("\t")
    if len(parts) < 5:
        return "no_data", proc.stdout.strip(), 0
    depth = int(parts[3])
    bases = parts[4]
    # For SNVs only — indels handled separately. Count alt-supporting bases.
    if len(ref) == 1 and len(alt) == 1:
        alt_upper = alt.upper()
        alt_count = bases.upper().count(alt_upper)
        # `.` and `,` are ref-matches (forward/reverse strand).
        ref_count = bases.count(".") + bases.count(",")
        if depth < min_depth:
            return "low_coverage", bases, depth
        if alt_count == 0:
            return "hom_ref", bases, depth
        # If the alt allele has support but DRAGEN didn't call it, flag it.
        return "unexpected_alt", bases, depth
    # For indels the pileup symbols are more complex; fall back to "low_coverage"
    # below threshold or "hom_ref" above.
    if depth < min_depth:
        return "low_coverage", bases, depth
    if "+" in bases or "-" in bases:
        return "unexpected_alt", bases, depth
    return "hom_ref", bases, depth


def lookup_in_vcf(sample_vf: pysam.VariantFile, sample_name: str,
                  contigs: set, chrom: str, pos: int, ref: str, alt: str):
    """Search the sample VCF for a record matching ref/alt at the position.

    Returns (gt_tuple, dp, sample_alt_matches) or (None, None, False).
    """
    chrom_candidates = []
    if chrom in contigs:
        chrom_candidates.append(chrom)
    alt_chrom = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
    if alt_chrom in contigs and alt_chrom not in chrom_candidates:
        chrom_candidates.append(alt_chrom)

    for c in chrom_candidates:
        try:
            it = sample_vf.fetch(c, pos - 1, pos)
        except (ValueError, KeyError):
            continue
        for rec in it:
            if rec.pos != pos or rec.ref != ref:
                continue
            for alt_in_rec in (rec.alts or []):
                if alt_in_rec == alt:
                    sample = (rec.samples.get(sample_name) or
                              list(rec.samples.values())[0])
                    gt = sample.get("GT")
                    dp = sample.get("DP")
                    return gt, dp, True
    return None, None, False


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--variants", type=Path, required=True, help="Variant YAML asset")
    ap.add_argument("--sample", type=Path, required=True, help="Indexed sample VCF")
    ap.add_argument("--cram", type=Path, default=None, help="Sample CRAM for pileup fallback")
    ap.add_argument("--reference-fasta", type=Path, default=None, help="Reference FASTA (needed with --cram)")
    ap.add_argument("--out", type=Path, required=True, help="Per-variant TSV")
    ap.add_argument("--out-combinations", type=Path, default=None,
                    help="Optional: per-combination TSV (APOE diplotype, HFE compound, etc.)")
    ap.add_argument("--sample-name", default=None)
    ap.add_argument("--min-depth", type=int, default=10,
                    help="Minimum CRAM depth to call hom_ref vs no_call (default 10)")
    args = ap.parse_args()

    if args.cram and not args.reference_fasta:
        ap.error("--cram requires --reference-fasta")

    with open(args.variants) as fh:
        asset = yaml.safe_load(fh)

    sample_vf = pysam.VariantFile(str(args.sample))
    if args.sample_name and args.sample_name in sample_vf.header.samples:
        sample_name = args.sample_name
    else:
        sample_name = list(sample_vf.header.samples)[0]
    contigs = set(sample_vf.header.contigs.keys())

    rows = []
    lookup_results: dict[str, dict] = {}  # keyed by variant id, used for combinations
    for v in asset["variants"]:
        gt, dp, found = lookup_in_vcf(
            sample_vf, sample_name, contigs,
            v["chrom"], v["pos"], v["ref"], v["alt"]
        )
        zyg = zygosity_class(gt) if found else None
        risk_count = None
        source = "vcf"
        pileup = ""
        if found:
            sample_alt_is_risk = (v["alt"] == v["risk_allele"])
            risk_count = risk_allele_count(zyg, gt, sample_alt_is_risk)
        elif args.cram:
            inference, pileup, depth = cram_pileup_inference(
                args.cram, args.reference_fasta,
                v["chrom"], v["pos"], v["ref"], v["alt"], args.min_depth,
            )
            source = f"cram:{inference}"
            dp = depth
            if inference == "hom_ref":
                zyg = "hom_ref"
                gt = (0, 0)
                sample_alt_is_risk = (v["alt"] == v["risk_allele"])
                risk_count = risk_allele_count(zyg, gt, sample_alt_is_risk)
            elif inference == "low_coverage":
                zyg = "low_coverage"
                risk_count = -1
            elif inference == "unexpected_alt":
                zyg = "unexpected_alt_in_pileup"
                risk_count = -1
            else:
                zyg = "no_data"
                risk_count = -1
        else:
            zyg = "no_vcf_record"
            risk_count = -1
        row = {
            "id": v["id"],
            "rsid": v["rsid"],
            "gene": v["gene"],
            "name": v["name"],
            "chrom": v["chrom"],
            "pos": v["pos"],
            "ref": v["ref"],
            "alt": v["alt"],
            "risk_allele": v["risk_allele"],
            "sample_gt": gt_summary(gt) if gt else "",
            "sample_zygosity": zyg,
            "sample_dp": dp if dp is not None else "",
            "risk_allele_count": risk_count,
            "source": source,
            "pileup": pileup,
            "category": v.get("category", ""),
            "evidence_tier": v.get("evidence_tier", ""),
            "condition": v.get("condition", ""),
            "interpretation_group": v.get("interpretation_group", ""),
            "interpretation": (v.get("interpretation") or "").strip().replace("\n", " "),
            "citations": "; ".join(v.get("citations", [])),
        }
        rows.append(row)
        lookup_results[v["id"]] = {"row": row, "raw": v}

    # Write per-variant TSV.
    if rows:
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    print(f"Wrote {len(rows)} variants → {args.out}", file=sys.stderr)

    # Handle combinations.
    if args.out_combinations and asset.get("combinations"):
        combo_rows = []
        for combo in asset["combinations"]:
            handler = combo.get("handler", "")
            handler_path = Path(__file__).parent / Path(handler).name
            if not handler_path.exists():
                combo_rows.append({
                    "id": combo["id"],
                    "description": combo["description"].strip().replace("\n", " "),
                    "interpretation": f"handler {handler_path} not found",
                    "components": ",".join(combo["component_variants"]),
                })
                continue
            # Pass component genotypes to the handler as JSON via stdin
            import json
            payload = {
                "components": {
                    cid: {
                        "rsid": lookup_results[cid]["row"]["rsid"],
                        "gene": lookup_results[cid]["row"]["gene"],
                        "ref": lookup_results[cid]["row"]["ref"],
                        "alt": lookup_results[cid]["row"]["alt"],
                        "risk_allele": lookup_results[cid]["row"]["risk_allele"],
                        "sample_gt": lookup_results[cid]["row"]["sample_gt"],
                        "sample_zygosity": lookup_results[cid]["row"]["sample_zygosity"],
                        "sample_dp": lookup_results[cid]["row"]["sample_dp"],
                    }
                    for cid in combo["component_variants"]
                    if cid in lookup_results
                }
            }
            proc = subprocess.run(
                [sys.executable, str(handler_path)],
                input=json.dumps(payload), capture_output=True, text=True, check=False,
            )
            if proc.returncode != 0:
                combo_rows.append({
                    "id": combo["id"],
                    "description": combo["description"].strip().replace("\n", " "),
                    "interpretation": f"handler error: {proc.stderr.strip()[:200]}",
                    "components": ",".join(combo["component_variants"]),
                })
                continue
            handler_out = json.loads(proc.stdout)
            combo_rows.append({
                "id": combo["id"],
                "description": combo["description"].strip().replace("\n", " "),
                "diplotype_or_classification": handler_out.get("classification", ""),
                "interpretation": handler_out.get("interpretation", ""),
                "components": ",".join(combo["component_variants"]),
                "component_summary": handler_out.get("component_summary", ""),
            })
        if combo_rows:
            with open(args.out_combinations, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(combo_rows[0].keys()), delimiter="\t")
                w.writeheader()
                w.writerows(combo_rows)
            print(f"Wrote {len(combo_rows)} combinations → {args.out_combinations}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
