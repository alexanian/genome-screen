process CURATED_LOOKUP {
    tag "${sample_vcf.baseName}"

    conda 'bioconda::pysam=0.22 bioconda::samtools=1.20 bioconda::htslib=1.20 conda-forge::pyyaml=6.0'

    input:
    path sample_vcf
    path sample_vcf_idx
    path sample_cram             // required for the coverage-fallback to distinguish hom_ref from no_call
    path sample_cram_idx
    path reference_fasta
    path reference_fasta_idx
    path variants_yaml           // either curated_risk_variants.yaml or nutrigenomic_variants.yaml

    output:
    path '04_curated_risk_variants.tsv',  emit: per_variant
    path '04_curated_combinations.tsv',   emit: combinations, optional: true

    script:
    """
    set -euo pipefail

    if [[ ! -s "${sample_cram}" || ! -s "${reference_fasta}" ]]; then
        echo "ERROR: curated lookup requires both a CRAM and a reference FASTA so that" >&2
        echo "       positions absent from the VCF can be distinguished as hom_ref vs no_call." >&2
        exit 1
    fi

    python3 ${projectDir}/bin/lookup_variants.py \\
        --variants "${variants_yaml}" \\
        --sample "${sample_vcf}" \\
        --cram "${sample_cram}" \\
        --reference-fasta "${reference_fasta}" \\
        --out 04_curated_risk_variants.tsv \\
        --out-combinations 04_curated_combinations.tsv \\
        --min-depth ${params.min_coverage_depth}
    """
}
