# genome-screen

A reproducible Nextflow pipeline that screens consumer-WGS data (Nebula, Dante, sequencing.com, etc.) for clinically actionable variants. Designed as a *pre-PRS* screen — polygenic scoring is handled separately by [pgsc_calc](https://github.com/PGScatalog/pgsc_calc).

## What it does

1. **Reference check** — confirms which genome build your VCF was called against (GRCh38, GRCh37, T2T).
2. **VCF QC** — ts/tv, het/hom, FILTER breakdown, sex inference.
3. **ACMG SF v3.2 secondary findings** — ClinVar Pathogenic/Likely-Pathogenic variants in the ACMG recommended genes.
4. **Other ClinVar P/LP** — same, genome-wide, ≥2-star review.
5. **Curated risk variants** — APOE ε2/ε3/ε4, F5 Leiden, F2, HFE, LRRK2 G2019S, etc. (variants ClinVar does not flag as pathogenic but have clear actionable evidence).
6. **Nutrigenomics** — variants graded ≥ B- by MyGeneFood (MTHFR, FUT2, LCT, ALDH2, FTO, TCF7L2, APOB, FADS1, …).
7. **Pharmacogenomics** — runs [PharmCAT](https://pharmcat.org/) for CPIC-aligned drug-gene recommendations.

All findings are consolidated into a single `report.html`.

## Quick start: reference check (no Nextflow required)

Before you commit to running the full pipeline, sanity-check that your VCF is the build you think it is:

```bash
python bin/check_reference.py /path/to/your_sample.vcf.gz
```

This single script has no Nextflow dependency and can be run against any VCF. It reads the header, fingerprints the contig MD5s, and tells you which reference panel was used (or warns if it doesn't recognize one). Useful as a triage step for any genome-analysis project.

## Full pipeline

Requires [Nextflow](https://www.nextflow.io/) ≥ 23.10 and either conda/mamba or Docker.

```bash
nextflow run . \
  --sample_vcf /path/to/sample.vcf.gz \
  --sample_cram /path/to/sample.cram \
  --reference_fasta /path/to/Homo_sapiens_assembly38.fasta \
  --expected_build GRCh38 \
  -profile conda    # or: -profile docker
```

Outputs are written to `results/`.

## Why these specific sources?

- **ClinVar ≥ 2-star** filters out single-submitter and conflicting calls, which are noisy for self-screening.
- **ACMG SF v3.2** is the consensus list of gene-disease pairs for which incidental findings are considered worth returning to patients (Miller DT et al. 2023, [PMID 37347242](https://pubmed.ncbi.nlm.nih.gov/37347242/)).
- **MyGeneFood science grades** (https://www.mygenefood.com/science-grade/) provide a published evidence rubric for nutrigenomic variants; we seed our list from grades A through B-, then manually verify cited PMIDs before promoting.
- **PharmCAT + CPIC** is the standard, peer-reviewed framework for translating genotype to drug recommendations.

## What it does *not* do

- **Polygenic risk scores** — use [pgsc_calc](https://github.com/PGScatalog/pgsc_calc) downstream. PRS requires re-calling at PGS Catalog positions (because effect alleles can equal the reference allele), which this pipeline deliberately skips. See [PGSC discussion #123](https://github.com/PGScatalog/pgsc_calc/discussions/123).
- **Structural variants / CNVs** — use Manta, Delly, or GATK gCNV.
- **Diagnostic interpretation** — this is a screening tool. Any positive finding should be confirmed by an accredited clinical lab and discussed with a genetic counselor.

## Privacy

Personal genomic data is *never* committed to this repo. Sample paths are passed in at runtime. The `.gitignore` excludes `*.vcf*`, `*.cram*`, and `results/` as a safety net.

## License

MIT for code. Curated assets cite their primary sources; see `assets/*.yaml` for per-entry citations.
