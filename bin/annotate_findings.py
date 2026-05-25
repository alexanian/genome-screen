#!/usr/bin/env python3
"""Cross-reference filtered ClinVar P/LP variants against a sample VCF and
emit two TSVs:

  - 02_acmg_sf_findings.tsv   — hits in ACMG SF v3.x genes (highest priority)
  - 03_clinvar_other_findings.tsv — P/LP hits in other genes

Each row carries the sample genotype, ClinVar annotation, and (for ACMG SF
hits) the per-gene reporting criterion + whether the sample's zygosity meets it.

Why split? ACMG SF genes have been judged actionable by an expert panel; other
P/LP hits are still informative but warrant more cautious interpretation.

Usage:
    python annotate_findings.py
        --sample sample.vcf.gz
        --clinvar clinvar.filtered.vcf.gz
        --acmg-sf assets/acmg_sf_v3_3.yaml
        --out-acmg-sf 02_acmg_sf_findings.tsv
        --out-other  03_clinvar_other_findings.tsv
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

try:
    import pysam
    import yaml
except ImportError as exc:
    print(f"ERROR: missing dependency ({exc}). Install via the genome-screen conda env.",
          file=sys.stderr)
    sys.exit(1)


# Sequence Ontology IDs for "truncating" molecular consequences (used by TTN
# and any other 'truncating_plp_only' genes).
TRUNCATING_SO_IDS = {
    "SO:0001587",  # stop_gained
    "SO:0001589",  # frameshift_variant
    "SO:0001574",  # splice_acceptor_variant
    "SO:0001575",  # splice_donor_variant
    "SO:0001572",  # exon_loss_variant
    "SO:0001578",  # stop_lost
    "SO:0001619",  # non_coding_transcript_variant — NOT truncating, listed for clarity (excluded)
}
TRUNCATING_SO_IDS.discard("SO:0001619")

# HFE c.845G>A / p.C282Y is rs1800562 on GRCh38 chr6:26092913 G>A.
HFE_C282Y = {"chrom": ("chr6", "6"), "pos": 26092913, "ref": "G", "alt": "A"}


def load_acmg(yaml_path: Path) -> dict:
    """Returns {gene_symbol: {inheritance, criteria, diseases}}."""
    with open(yaml_path) as fh:
        d = yaml.safe_load(fh)
    out = {}
    for g in d["genes"]:
        out[g["symbol"]] = {
            "inheritance": g["inheritance"],
            "added_in_version": g.get("added_in_version"),
            "criteria": sorted({dis["reporting_criterion"] for dis in g["diseases"]}),
            "diseases": [dis["name"] for dis in g["diseases"]],
        }
    return out


def parse_geneinfo(info_value: str | None) -> list[str]:
    """ClinVar GENEINFO is 'GENE:ID|GENE2:ID2'. Returns symbol list, first first."""
    if not info_value:
        return []
    return [pair.split(":", 1)[0] for pair in info_value.split("|") if pair]


def is_truncating(mc_value: str | None) -> bool:
    """ClinVar MC is 'SO:NNNNN|consequence,SO:NNNNN|consequence'."""
    if not mc_value:
        return False
    for entry in mc_value.split(","):
        so_id = entry.split("|", 1)[0]
        if so_id in TRUNCATING_SO_IDS:
            return True
    return False


def gt_summary(gt_tuple) -> str:
    """Render a pysam genotype tuple as a readable string ('0/1', '1/1', './.')."""
    if gt_tuple is None or any(a is None for a in gt_tuple):
        return "./."
    return "/".join(str(a) for a in gt_tuple)


def zygosity_class(gt_tuple) -> str:
    """Classify GT into hom_ref / het / hom_alt / hemi / nocall."""
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


def lookup_sample_genotype(sample_vcf: pysam.VariantFile, sample_name: str,
                            chrom_candidates: list[str], pos: int,
                            ref: str, alt: str):
    """Search the sample VCF for a record matching ref/alt at the given position.

    Returns (record, gt_tuple, dp) or (None, None, None) if no match.
    """
    for chrom in chrom_candidates:
        try:
            it = sample_vcf.fetch(chrom, pos - 1, pos)
        except (ValueError, KeyError):
            continue
        for rec in it:
            if rec.pos != pos:
                continue
            if rec.ref != ref:
                continue
            # ALT may be multiallelic in the sample
            for i, sample_alt in enumerate(rec.alts or [], start=1):
                if sample_alt == alt:
                    sample = rec.samples.get(sample_name) or list(rec.samples.values())[0]
                    gt = sample.get("GT")
                    # Translate multiallelic alt index into the GT call: pysam already
                    # gives the correct allele index, so we just return the genotype.
                    dp = sample.get("DP")
                    return rec, gt, dp
    return None, None, None


def evaluate_criterion(criterion: str, zyg: str, hits_in_gene_this_sample: int,
                       *, is_truncating_v: bool, is_hfe_c282y: bool) -> tuple[bool, str]:
    """Decide whether the sample's zygosity meets the ACMG SF reporting criterion.

    Returns (meets, rationale).
    """
    if criterion == "all_plp":
        if zyg in ("het", "hom_alt", "hemi_alt"):
            return True, f"any P/LP carrier ({zyg})"
        return False, f"variant not present in sample ({zyg})"
    if criterion == "any_zygosity_plp":
        if zyg in ("het", "hom_alt", "hemi_alt"):
            return True, f"any zygosity actionable ({zyg})"
        return False, f"variant not present ({zyg})"
    if criterion == "hemi_homo_or_2het_plp":
        if zyg in ("hom_alt", "hemi_alt"):
            return True, f"homozygous / hemizygous ({zyg})"
        if zyg == "het" and hits_in_gene_this_sample >= 2:
            return True, "compound het: ≥2 P/LP variants in this gene"
        return False, f"single het not actionable here ({zyg})"
    if criterion == "two_plp_variants":
        if zyg == "hom_alt":
            return True, "homozygous P/LP"
        if zyg == "het" and hits_in_gene_this_sample >= 2:
            return True, "compound het: ≥2 P/LP variants in this gene"
        return False, "needs ≥2 P/LP variants"
    if criterion == "truncating_plp_only":
        if not is_truncating_v:
            return False, "ClinVar consequence not truncating"
        if zyg in ("het", "hom_alt", "hemi_alt"):
            return True, f"truncating variant present ({zyg})"
        return False, "truncating variant absent"
    if criterion == "hfe_c282y_homozygotes":
        if not is_hfe_c282y:
            return False, "not the p.C282Y variant"
        if zyg == "hom_alt":
            return True, "p.C282Y homozygous"
        return False, f"only homozygotes reported ({zyg})"
    return False, f"unknown criterion {criterion!r}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True, type=Path, help="Indexed sample VCF (vcf.gz + tbi)")
    ap.add_argument("--clinvar", required=True, type=Path, help="Filtered ClinVar VCF (vcf.gz + tbi)")
    ap.add_argument("--acmg-sf", required=True, type=Path, help="ACMG SF YAML asset")
    ap.add_argument("--out-acmg-sf", required=True, type=Path)
    ap.add_argument("--out-other", required=True, type=Path)
    ap.add_argument("--sample-name", default=None,
                    help="Sample name in VCF (defaults to the first sample)")
    args = ap.parse_args()

    acmg = load_acmg(args.acmg_sf)
    sample_vf = pysam.VariantFile(str(args.sample))
    if args.sample_name and args.sample_name in sample_vf.header.samples:
        sample_name = args.sample_name
    else:
        sample_name = list(sample_vf.header.samples)[0]
        print(f"Using sample '{sample_name}' (first in VCF).", file=sys.stderr)

    sample_contigs = set(sample_vf.header.contigs.keys())

    # Pass 1: walk ClinVar, look up matching sample variants, collect raw hits.
    raw_hits = []
    clinvar_vf = pysam.VariantFile(str(args.clinvar))
    n_cv = 0
    for rec in clinvar_vf:
        n_cv += 1
        if n_cv % 10000 == 0:
            print(f"  ...scanned {n_cv:,} ClinVar records", file=sys.stderr)
        ref = rec.ref
        for alt in (rec.alts or []):
            # Build candidate chromosome names for the sample VCF.
            chrom_candidates = []
            if rec.chrom in sample_contigs:
                chrom_candidates.append(rec.chrom)
            alt_chrom = rec.chrom[3:] if rec.chrom.startswith("chr") else f"chr{rec.chrom}"
            if alt_chrom in sample_contigs and alt_chrom not in chrom_candidates:
                chrom_candidates.append(alt_chrom)
            if not chrom_candidates:
                continue
            _, gt, dp = lookup_sample_genotype(
                sample_vf, sample_name, chrom_candidates, rec.pos, ref, alt
            )
            if gt is None:
                continue
            zyg = zygosity_class(gt)
            if zyg in ("hom_ref", "hemi_ref", "nocall"):
                continue
            genes = parse_geneinfo(rec.info.get("GENEINFO"))
            primary_gene = genes[0] if genes else ""
            raw_hits.append({
                "chrom": rec.chrom,
                "pos": rec.pos,
                "ref": ref,
                "alt": alt,
                "rsid": rec.id or "",
                "sample_gt": gt_summary(gt),
                "sample_zygosity": zyg,
                "sample_dp": dp if dp is not None else "",
                "gene": primary_gene,
                "all_genes": "|".join(genes),
                "clnsig": rec.info.get("CLNSIG", ("",))[0] if isinstance(rec.info.get("CLNSIG"), tuple) else rec.info.get("CLNSIG", ""),
                "clnrevstat": rec.info.get("CLNREVSTAT", ("",))[0] if isinstance(rec.info.get("CLNREVSTAT"), tuple) else rec.info.get("CLNREVSTAT", ""),
                "condition": ";".join(rec.info.get("CLNDN", []) or []) if rec.info.get("CLNDN") else "",
                "allele_id": rec.info.get("ALLELEID", ""),
                "clnhgvs": ";".join(rec.info.get("CLNHGVS", []) or []) if rec.info.get("CLNHGVS") else "",
                "mc": ";".join(rec.info.get("MC", []) or []) if rec.info.get("MC") else "",
                "clnvc": rec.info.get("CLNVC", ""),
            })

    print(f"ClinVar records scanned: {n_cv:,}", file=sys.stderr)
    print(f"Sample hits found:       {len(raw_hits):,}", file=sys.stderr)

    # Pass 2: count hits per ACMG SF gene for compound-het logic.
    hits_per_gene = defaultdict(int)
    for h in raw_hits:
        if h["gene"] in acmg:
            hits_per_gene[h["gene"]] += 1

    # Pass 3: classify each hit, apply ACMG SF criterion.
    acmg_rows = []
    other_rows = []
    for h in raw_hits:
        if h["gene"] in acmg:
            gene_info = acmg[h["gene"]]
            criteria = gene_info["criteria"]
            # If a gene has multiple criteria (e.g. SCN5A is all_plp for 3 diseases),
            # apply the most permissive one for the meets-criterion flag, but list all.
            n_in_gene = hits_per_gene[h["gene"]]
            results = []
            for crit in criteria:
                ok, why = evaluate_criterion(
                    crit, h["sample_zygosity"], n_in_gene,
                    is_truncating_v=is_truncating(h["mc"]),
                    is_hfe_c282y=(
                        h["chrom"] in HFE_C282Y["chrom"]
                        and h["pos"] == HFE_C282Y["pos"]
                        and h["ref"] == HFE_C282Y["ref"]
                        and h["alt"] == HFE_C282Y["alt"]
                    ),
                )
                results.append((crit, ok, why))
            meets = any(ok for _, ok, _ in results)
            crit_summary = ";".join(f"{c}={'Y' if ok else 'N'}({why})" for c, ok, why in results)
            acmg_rows.append({
                **h,
                "acmg_inheritance": gene_info["inheritance"],
                "acmg_added_in": gene_info["added_in_version"],
                "acmg_criteria_evaluated": ",".join(criteria),
                "meets_criterion": "Y" if meets else "N",
                "criterion_detail": crit_summary,
            })
        else:
            other_rows.append(h)

    # Write outputs.
    def write_tsv(rows: list[dict], path: Path):
        if not rows:
            # Still write the header so downstream parsers don't choke
            path.write_text("# no findings\n")
            return
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    write_tsv(acmg_rows, args.out_acmg_sf)
    write_tsv(other_rows, args.out_other)

    print(f"ACMG SF findings:  {len(acmg_rows):,}  → {args.out_acmg_sf}", file=sys.stderr)
    print(f"Other P/LP hits:   {len(other_rows):,}  → {args.out_other}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
