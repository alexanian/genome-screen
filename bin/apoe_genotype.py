#!/usr/bin/env python3
"""Derive APOE ε2/ε3/ε4 diplotype from rs429358 + rs7412 genotypes.

Called by lookup_variants.py via JSON-on-stdin / JSON-on-stdout. The two
component variants are passed in as a {"components": {...}} dict.

The diplotype-derivation table is the standard textbook one:

  rs429358 GT  | rs7412 GT  | Diplotype
  -------------+------------+----------
  TT (ref/ref) | CC (ref/ref) | ε3/ε3 (most common)
  TC (het)     | CC (ref/ref) | ε3/ε4
  CC (hom-alt) | CC (ref/ref) | ε4/ε4
  TT (ref/ref) | CT (het)     | ε2/ε3
  TC (het)     | CT (het)     | ε2/ε4 (ε1 alleles essentially extinct)
  TT (ref/ref) | TT (hom-alt) | ε2/ε2
  others       | others       | impossible (ε1 carrier) — flagged as inconsistent

Note: this assumes the C allele at rs429358 and the T allele at rs7412 are
in trans when both are heterozygous (ε2/ε4 rather than ε1/ε3). The ε1
allele is so rare (<0.01% AF in all studied populations) that this is
effectively always correct.
"""

from __future__ import annotations

import json
import sys

DIPLOTYPE_TABLE = {
    ("hom_ref", "hom_ref"): ("ε3/ε3", "Most common; no APOE-mediated risk."),
    ("het",     "hom_ref"): ("ε3/ε4", "Heterozygous ε4: ~2-3x AD risk, ~+10% LDL."),
    ("hom_alt", "hom_ref"): ("ε4/ε4", "Homozygous ε4: ~10-12x AD risk, elevated CAD risk."),
    ("hom_ref", "het"):     ("ε2/ε3", "Heterozygous ε2: mildly AD-protective."),
    ("het",     "het"):     ("ε2/ε4", "Both ε2 and ε4 — mixed effect; ε4 risk usually dominates for AD."),
    ("hom_ref", "hom_alt"): ("ε2/ε2", "Homozygous ε2: AD-protective; necessary background for type-III hyperlipoproteinemia (penetrance ~10%)."),
}


def classify(zyg_429358: str, zyg_7412: str) -> tuple[str, str]:
    key = (zyg_429358, zyg_7412)
    if key in DIPLOTYPE_TABLE:
        return DIPLOTYPE_TABLE[key]
    if "nocall" in key or "no_data" in key or "low_coverage" in key:
        return ("indeterminate",
                f"Cannot determine APOE diplotype — rs429358={zyg_429358}, rs7412={zyg_7412}. "
                "Confirm coverage at both positions.")
    return ("inconsistent",
            f"Genotype combination {zyg_429358}/{zyg_7412} would imply a rare ε1 allele "
            "or a genotyping artifact. Confirm both calls manually.")


def main() -> int:
    payload = json.load(sys.stdin)
    comps = payload.get("components", {})
    rs429358 = comps.get("APOE_rs429358", {})
    rs7412 = comps.get("APOE_rs7412", {})
    if not rs429358 or not rs7412:
        print(json.dumps({
            "classification": "incomplete",
            "interpretation": "Missing one or both APOE component genotypes.",
            "component_summary": "",
        }))
        return 0

    diplotype, interp = classify(rs429358["sample_zygosity"], rs7412["sample_zygosity"])
    component_summary = (
        f"rs429358={rs429358['sample_gt']} ({rs429358['sample_zygosity']}, "
        f"DP={rs429358.get('sample_dp') or '?'}); "
        f"rs7412={rs7412['sample_gt']} ({rs7412['sample_zygosity']}, "
        f"DP={rs7412.get('sample_dp') or '?'})"
    )
    print(json.dumps({
        "classification": diplotype,
        "interpretation": interp,
        "component_summary": component_summary,
    }))
    return 0


if __name__ == "__main__":
    sys.exit(main())
