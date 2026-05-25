// Nutrigenomic variant lookup — reuses the same lookup_variants.py
// machinery as the curated-risk-variants step, just with a different
// YAML and output filename. The script handles per-variant VCF lookup
// with CRAM-pileup fallback for positions absent from the sample VCF.

process NUTRIGENOMICS {
    tag "${sample_vcf.baseName}"

    conda 'bioconda::pysam=0.22 bioconda::samtools=1.20 bioconda::htslib=1.20 conda-forge::pyyaml=6.0'

    input:
    path sample_vcf
    path sample_vcf_idx
    path sample_cram
    path sample_cram_idx
    path reference_fasta
    path reference_fasta_idx
    path nutrigenomic_yaml

    output:
    path '05_nutrigenomics.tsv', emit: per_variant

    script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/lookup_variants.py \\
        --variants "${nutrigenomic_yaml}" \\
        --sample "${sample_vcf}" \\
        --cram "${sample_cram}" \\
        --reference-fasta "${reference_fasta}" \\
        --out 05_nutrigenomics.tsv \\
        --min-depth ${params.min_coverage_depth}
    """
}
