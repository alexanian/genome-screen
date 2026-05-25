#!/usr/bin/env python3
"""Classify HFE hemochromatosis status from joint C282Y + H63D genotypes.

ACMG SF reports ONLY C282Y/C282Y homozygotes. This handler covers the rest
of the clinically-relevant combinations the SF list deliberately excludes:

  C282Y zygosity | H63D zygosity | Class                  | Risk
  ---------------+---------------+------------------------+-------------------
  hom_alt        | hom_ref       | C282Y/C282Y            | High (severe HH)
  het            | het           | C282Y/H63D compound    | Moderate (lower
                                  |  het (assumed trans)   |  penetrance HH)
  het            | hom_ref       | C282Y carrier          | Mild
  hom_ref        | hom_alt       | H63D/H63D              | Low (mild iron
                                                              load possible)
  hom_ref        | het           | H63D carrier           | Negligible
  hom_ref        | hom_ref       | Wildtype               | None
  hom_alt        | het/hom_alt   | C282Y hom + H63D       | Rare; treat as
                                                              C282Y/C282Y for
                                                              risk purposes

Assumption: het/het is reported as compound heterozygous (in trans), which
is the much more common configuration. The in-cis case (both variants on
the same chromosome) is rare but cannot be determined without phasing;
flag this in the interpretation for the user to consider.
"""

from __future__ import annotations

import json
import sys


def classify(c282y: str, h63d: str) -> tuple[str, str]:
    if c282y == "hom_alt":
        if h63d == "hom_ref":
            return ("C282Y/C282Y homozygous",
                    "Highest-penetrance hereditary hemochromatosis (HH) genotype. "
                    "This is what the ACMG SF list reports. Confirm with serum ferritin + "
                    "transferrin saturation; treatment is straightforward (phlebotomy).")
        return ("C282Y hom + H63D variant",
                "Atypical: C282Y homozygous with additional H63D allele. Treat clinically "
                "as C282Y/C282Y for HH risk purposes.")
    if c282y == "het":
        if h63d == "het":
            return ("C282Y/H63D compound heterozygous (assumed)",
                    "Compound heterozygous configuration if the alleles are on different "
                    "chromosomes (in trans), which is the common case. Moderate iron-overload "
                    "risk, lower penetrance than C282Y/C282Y. Reasonable to check serum "
                    "ferritin + transferrin saturation every 5-10 years. Cannot distinguish "
                    "from in-cis (both on same chromosome) without parental data or phasing.")
        if h63d == "hom_alt":
            return ("C282Y/H63D + H63D",
                    "C282Y heterozygote plus H63D homozygote. Effectively similar to "
                    "compound het for iron-overload risk; check baseline iron studies.")
        return ("C282Y carrier (heterozygous)",
                "Single C282Y allele. Generally asymptomatic; minor iron-handling "
                "differences possible. Carrier status relevant for offspring if partner is "
                "also a carrier.")
    if c282y == "hom_ref":
        if h63d == "hom_alt":
            return ("H63D/H63D homozygous",
                    "Mild — most carriers do not develop clinical hemochromatosis without "
                    "additional cofactors (heavy alcohol, hepatitis C). No routine "
                    "surveillance indicated unless iron studies are abnormal.")
        if h63d == "het":
            return ("H63D carrier (heterozygous)",
                    "Negligible HH risk in isolation.")
        return ("Wildtype",
                "No hemochromatosis-associated HFE variants detected.")
    if "nocall" in (c282y, h63d) or "low_coverage" in (c282y, h63d) or "no_data" in (c282y, h63d):
        return ("indeterminate",
                f"Cannot classify — C282Y={c282y}, H63D={h63d}. Confirm coverage.")
    return ("indeterminate", f"Unrecognized state C282Y={c282y}, H63D={h63d}.")


def main() -> int:
    payload = json.load(sys.stdin)
    comps = payload.get("components", {})
    c282y = comps.get("HFE_C282Y", {})
    h63d = comps.get("HFE_H63D", {})
    if not c282y or not h63d:
        print(json.dumps({
            "classification": "incomplete",
            "interpretation": "Missing one or both HFE component genotypes.",
            "component_summary": "",
        }))
        return 0
    classification, interp = classify(c282y["sample_zygosity"], h63d["sample_zygosity"])
    component_summary = (
        f"C282Y={c282y['sample_gt']} ({c282y['sample_zygosity']}, "
        f"DP={c282y.get('sample_dp') or '?'}); "
        f"H63D={h63d['sample_gt']} ({h63d['sample_zygosity']}, "
        f"DP={h63d.get('sample_dp') or '?'})"
    )
    print(json.dumps({
        "classification": classification,
        "interpretation": interp,
        "component_summary": component_summary,
    }))
    return 0


if __name__ == "__main__":
    sys.exit(main())
