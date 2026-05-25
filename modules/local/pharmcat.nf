// PharmCAT — CPIC-graded pharmacogenomic recommendations from a sample VCF.
//
// Uses the official PharmCAT pipeline (Python wrapper + Java JAR + position VCF)
// rather than running PharmCAT and the preprocessor as separate processes — the
// bundled pipeline handles preprocessing, allele matching, phenotyping, and
// report generation in one step.
//
// Inputs:
//   - sample_vcf: variant-sites VCF (passed --absent-to-ref so PharmCAT treats
//     PGx positions absent from the VCF as homozygous reference; required for
//     non-gVCF input like Nebula's DRAGEN-called VCFs)
//
// Outputs:
//   - <sample>.report.html       — human-readable PharmCAT report
//   - <sample>.phenotype.json    — gene-level diplotype + phenotype JSON
//   - <sample>.match.json        — Named Allele Matcher detail
//   - <sample>.preprocessed.vcf.bgz — preprocessed VCF restricted to PGx positions
//
// Known limitations (short-read WGS):
//   - CYP2D6 typically reported as "No Result" — needs Cyrius/Stargazer for
//     structural-variant-aware calling. The single largest PGx gap.
//   - HLA-A / HLA-B "No Result" — needs HLA-LA or Optitype.
//   - MT-RNR1 may not be called if mtDNA coverage is unreliable.

process PHARMCAT {
    tag "${sample_vcf.simpleName}"

    conda 'bioconda::bcftools=1.20 bioconda::htslib=1.20 conda-forge::openjdk=17 conda-forge::pandas conda-forge::colorama conda-forge::packaging'

    input:
    path sample_vcf
    path sample_vcf_idx
    path pharmcat_dir           // unpacked pharmcat-pipeline release directory

    output:
    path "*.report.html",         emit: report
    path "*.phenotype.json",      emit: phenotype
    path "*.match.json",          emit: match
    path "*.preprocessed.vcf.bgz",emit: preprocessed_vcf, optional: true
    path "*.missing_pgx_var.vcf", emit: missing_pgx,      optional: true

    script:
    """
    set -euo pipefail
    python3 ${pharmcat_dir}/pharmcat_pipeline \\
        ${sample_vcf} \\
        --absent-to-ref \\
        -o .
    """
}
