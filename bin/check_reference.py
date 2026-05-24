#!/usr/bin/env python3
"""Identify which reference genome build a VCF was called against.

Reads only the VCF header — does not require the FASTA, the index, or bcftools.
Designed to be runnable standalone (no pipeline dependencies) so anyone can
sanity-check their consumer-WGS file before plugging it into a downstream tool.

Usage:
    python check_reference.py sample.vcf.gz [--signatures path/to/reference_signatures.yaml]
                                            [--expected GRCh38]
                                            [--json]

Exits 0 if a single best match is found and (if --expected is given) matches it.
Exits 2 if no signature matches confidently or the build differs from --expected.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import re
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    print(
        "ERROR: PyYAML is required. Install with `pip install pyyaml`.",
        file=sys.stderr,
    )
    sys.exit(1)


DEFAULT_SIGNATURES = Path(__file__).resolve().parent.parent / "assets" / "reference_signatures.yaml"

# Confidence thresholds (out of 100).
# A perfect chromosome-length match across the 5 primaries is worth 50 points;
# naming convention is 10; contig-count matches add up to 30; MD5s (when
# present) add 10 each. So a typical strong match scores ~80-100.
MIN_CONFIDENT_SCORE = 50
MIN_MATCH_DELTA = 15  # best must beat runner-up by this margin


# --- VCF header parsing ---------------------------------------------------

_CONTIG_RE = re.compile(r"^##contig=<([^>]+)>")
_REFERENCE_RE = re.compile(r"^##reference=(.+)$")


def open_vcf(path: Path) -> io.TextIOBase:
    """Open a VCF, gzipped or plain, as a text stream."""
    if str(path).endswith(".gz") or str(path).endswith(".bgz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_vcf_header(path: Path) -> dict:
    """Extract the bits of the VCF header relevant to reference identification.

    Returns: {
        'reference_line': str | None,
        'contigs': [{'name': str, 'length': int|None, 'md5': str|None}, ...],
    }
    """
    contigs = []
    reference_line = None
    with open_vcf(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.startswith("#"):
                break
            if line.startswith("#CHROM"):
                break
            ref_match = _REFERENCE_RE.match(line)
            if ref_match:
                reference_line = ref_match.group(1).strip()
                continue
            contig_match = _CONTIG_RE.match(line)
            if not contig_match:
                continue
            attrs = {}
            # Naive parse: split on commas not inside quotes/angle brackets.
            for kv in contig_match.group(1).split(","):
                if "=" not in kv:
                    continue
                k, _, v = kv.partition("=")
                attrs[k.strip()] = v.strip().strip('"')
            try:
                length = int(attrs.get("length", "")) if attrs.get("length") else None
            except ValueError:
                length = None
            contigs.append(
                {
                    "name": attrs.get("ID", ""),
                    "length": length,
                    "md5": attrs.get("M5") or attrs.get("md5"),
                }
            )
    return {"reference_line": reference_line, "contigs": contigs}


def derive_features(header: dict) -> dict:
    """Compute the features used to score against the signatures DB."""
    contigs = header["contigs"]
    by_name = {c["name"]: c for c in contigs}

    chr_prefix = sum(1 for c in contigs if c["name"].startswith("chr"))
    ensembl_naming = sum(1 for c in contigs if c["name"] in {"1", "2", "MT"} or c["name"].isdigit())
    naming = "chr_prefixed" if chr_prefix > ensembl_naming else "ensembl"

    # Count only true "_alt" alternate-locus contigs. HLA and decoy contigs are
    # categorised separately because the Broad bundle ships ~525 HLA contigs
    # which would otherwise dwarf the 261-_alt count and break matching.
    alt_count = sum(1 for c in contigs if c["name"].endswith("_alt"))

    return {
        "total_contigs": len(contigs),
        "alt_contigs": alt_count,
        "naming": naming,
        "by_name": by_name,
        "reference_line": header.get("reference_line"),
    }


# --- Scoring --------------------------------------------------------------


def score_signature(features: dict, sig: dict) -> tuple[int, list[str]]:
    """Score how well the parsed VCF matches one reference signature.

    Returns (score, notes). Higher is better; max is ~110 with all MD5s.
    """
    score = 0
    notes: list[str] = []

    # Naming convention (10 pts)
    if features["naming"] == sig.get("naming"):
        score += 10
        notes.append(f"naming={sig['naming']} ✓")
    else:
        notes.append(f"naming mismatch (vcf={features['naming']}, sig={sig.get('naming')})")

    # Primary-chromosome length match (up to 50 pts, 10 each)
    primaries = sig.get("primary_chromosomes", {})
    matched_primaries = 0
    for chrom, expected in primaries.items():
        observed = features["by_name"].get(chrom)
        if observed is None:
            notes.append(f"{chrom} absent from VCF header")
            continue
        if observed["length"] == expected["length"]:
            matched_primaries += 1
            score += 10
            # MD5 bonus (10 pts each) when both sides have one
            if expected.get("md5") and observed.get("md5") and expected["md5"] == observed["md5"]:
                score += 10
                notes.append(f"{chrom} length+MD5 ✓")
        else:
            notes.append(
                f"{chrom} length mismatch (vcf={observed['length']}, sig={expected['length']})"
            )
    if matched_primaries == len(primaries):
        notes.append(f"all {matched_primaries} primary chromosomes match")

    # Total contig count (15 pts, partial credit by ratio)
    sig_total = sig.get("total_contigs")
    if sig_total:
        if features["total_contigs"] == sig_total:
            score += 15
            notes.append(f"total_contigs={sig_total} ✓")
        else:
            # Don't penalize — VCFs sometimes drop contigs with no variants,
            # but matching exactly is a strong positive signal.
            notes.append(
                f"total_contigs differ (vcf={features['total_contigs']}, sig={sig_total})"
            )

    # Alt count (15 pts)
    sig_alts = sig.get("alt_contigs")
    if sig_alts is not None:
        if features["alt_contigs"] == sig_alts:
            score += 15
            notes.append(f"alt_contigs={sig_alts} ✓")
        elif sig_alts == 0 and features["alt_contigs"] > 0:
            notes.append(f"VCF has {features['alt_contigs']} alts; sig expects 0")
        # Otherwise: silent partial mismatch.

    # Extra-contig tiebreaker (10 pts each): some signatures declare
    # diagnostic decoy contigs (e.g. hs37d5) that uniquely identify them.
    for extra in sig.get("extra_contigs", []):
        if extra in features["by_name"]:
            score += 10
            notes.append(f"diagnostic contig {extra!r} present ✓")

    # Reference-line keyword tiebreaker (5 pts each): a #reference= URL
    # mentioning the signature variant is weak but useful corroboration.
    ref_line = (features.get("reference_line") or "").lower()
    if ref_line:
        for keyword in [sig.get("id"), sig.get("variant", "").split()[0].lower()]:
            if keyword and keyword.lower() in ref_line:
                score += 5
                notes.append(f"##reference mentions {keyword!r}")
                break

    return score, notes


# --- Main -----------------------------------------------------------------


def identify(vcf_path: Path, sig_path: Path) -> dict:
    header = parse_vcf_header(vcf_path)
    features = derive_features(header)
    with open(sig_path) as fh:
        signatures = yaml.safe_load(fh)["references"]

    scored = []
    for sig in signatures:
        score, notes = score_signature(features, sig)
        scored.append({"sig": sig, "score": score, "notes": notes})
    scored.sort(key=lambda r: r["score"], reverse=True)

    best = scored[0]
    runner_up = scored[1] if len(scored) > 1 else None

    # Two levels of confidence:
    #   - variant_confident: we know the exact sub-variant (e.g. hs37d5 vs b37)
    #   - build_confident:   we know the major build (GRCh38 vs GRCh37 vs T2T)
    # The build call is the one users typically care about; sub-variant ties
    # between e.g. b37 and hs37d5 should not cause an overall FAIL.
    variant_confident = best["score"] >= MIN_CONFIDENT_SCORE and (
        runner_up is None or (best["score"] - runner_up["score"]) >= MIN_MATCH_DELTA
    )
    next_diff_build = next(
        (r for r in scored[1:] if r["sig"]["build"] != best["sig"]["build"]),
        None,
    )
    build_confident = best["score"] >= MIN_CONFIDENT_SCORE and (
        next_diff_build is None
        or (best["score"] - next_diff_build["score"]) >= MIN_MATCH_DELTA
    )

    return {
        "vcf": str(vcf_path),
        "observed": {
            "total_contigs": features["total_contigs"],
            "alt_contigs": features["alt_contigs"],
            "naming": features["naming"],
            "reference_line": features["reference_line"],
        },
        "best_match": {
            "id": best["sig"]["id"],
            "build": best["sig"]["build"],
            "variant": best["sig"]["variant"],
            "score": best["score"],
            "notes": best["notes"],
        },
        "runner_up": {
            "id": runner_up["sig"]["id"],
            "score": runner_up["score"],
        }
        if runner_up
        else None,
        "variant_confident": variant_confident,
        "build_confident": build_confident,
        "confident": variant_confident,  # back-compat
        "all_scores": [{"id": r["sig"]["id"], "score": r["score"]} for r in scored],
    }


def format_report(result: dict, expected: str | None) -> str:
    obs = result["observed"]
    best = result["best_match"]
    lines = []
    lines.append(f"VCF: {result['vcf']}")
    lines.append("")
    lines.append("Observed header features:")
    lines.append(f"  reference= line   : {obs['reference_line'] or '(none)'}")
    lines.append(f"  naming convention : {obs['naming']}")
    lines.append(f"  total contigs     : {obs['total_contigs']}")
    lines.append(f"  alt contigs       : {obs['alt_contigs']}")
    lines.append("")
    if result["variant_confident"]:
        lines.append(f"RESULT: {best['build']} — {best['variant']}  (score={best['score']})")
    elif result["build_confident"]:
        lines.append(
            f"RESULT: {best['build']} (sub-variant uncertain — best guess {best['variant']}, score={best['score']})"
        )
        if result["runner_up"]:
            ru = result["runner_up"]
            lines.append(f"        Runner-up: {ru['id']} (score={ru['score']})")
    else:
        lines.append(
            f"RESULT: UNCERTAIN. Best guess: {best['build']} — {best['variant']} (score={best['score']})"
        )
        if result["runner_up"]:
            ru = result["runner_up"]
            lines.append(f"        Runner-up: {ru['id']} (score={ru['score']})")
    lines.append("")
    lines.append("Match notes:")
    for n in best["notes"]:
        lines.append(f"  - {n}")
    lines.append("")
    lines.append("All signatures scored:")
    for s in result["all_scores"]:
        lines.append(f"  {s['id']:20s}  {s['score']}")

    if expected:
        actual_build = best["build"]
        if actual_build == expected and result["build_confident"]:
            lines.append("")
            lines.append(f"EXPECTED-BUILD CHECK: PASS  (expected={expected}, observed={actual_build})")
        else:
            lines.append("")
            lines.append(
                f"EXPECTED-BUILD CHECK: FAIL  (expected={expected}, observed={actual_build}, build_confident={result['build_confident']})"
            )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Identify reference genome build from a VCF header."
    )
    parser.add_argument("vcf", type=Path, help="Path to a VCF or VCF.gz file")
    parser.add_argument(
        "--signatures",
        type=Path,
        default=DEFAULT_SIGNATURES,
        help=f"Path to reference_signatures.yaml (default: {DEFAULT_SIGNATURES})",
    )
    parser.add_argument(
        "--expected",
        type=str,
        default=None,
        help="Expected build (e.g. GRCh38). Exits 2 if mismatch.",
    )
    parser.add_argument("--json", action="store_true", help="Emit JSON instead of text")
    args = parser.parse_args()

    if not args.vcf.exists():
        print(f"ERROR: VCF not found: {args.vcf}", file=sys.stderr)
        return 1
    if not args.signatures.exists():
        print(f"ERROR: signatures file not found: {args.signatures}", file=sys.stderr)
        return 1

    result = identify(args.vcf, args.signatures)

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print(format_report(result, args.expected))

    if args.expected:
        if result["best_match"]["build"] != args.expected or not result["build_confident"]:
            return 2
    elif not result["build_confident"]:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
