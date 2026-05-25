process CLINVAR_INTERSECT {
    tag "${sample_vcf.baseName}"

    conda 'bioconda::pysam=0.22 bioconda::htslib=1.20 conda-forge::pyyaml=6.0'

    input:
    path sample_vcf
    path sample_vcf_idx
    path sample_cram             // optional; pass an empty file if absent
    path sample_cram_idx
    path reference_fasta         // required when CRAM provided (CRAM needs ref to decode)
    path reference_fasta_idx
    tuple path(clinvar_vcf), path(clinvar_tbi), path(clinvar_md5), path(clinvar_date)
    path acmg_sf_yaml

    output:
    path '02_acmg_sf_findings.tsv',                  emit: acmg_sf
    path '02_acmg_sf_findings_verification.sh',      emit: acmg_sf_verify,  optional: true
    path '03_clinvar_other_findings.tsv',            emit: other
    path '03_clinvar_other_findings_verification.sh',emit: other_verify,    optional: true
    path 'clinvar.filtered.vcf.gz',                  emit: filtered_vcf
    path 'clinvar.filtered.vcf.gz.tbi',              emit: filtered_tbi

    script:
    // ClinVar's GRCh38 VCF uses unprefixed contigs ('1' not 'chr1'). Detect
    // the sample's naming convention and rewrite ClinVar contigs if needed.
    """
    set -euo pipefail

    # Decide whether to add the chr prefix when filtering ClinVar.
    SAMPLE_FIRST_CONTIG=\$(zcat "${sample_vcf}" | grep -m1 "^##contig=<ID=" | sed -E 's/^##contig=<ID=([^,>]+).*/\\1/')
    if [[ "\${SAMPLE_FIRST_CONTIG}" == chr* ]]; then
        RENAME_FLAG="--rename-chrs"
    else
        RENAME_FLAG=""
    fi
    echo "Sample uses contig naming: \${SAMPLE_FIRST_CONTIG} → rename flag: \${RENAME_FLAG:-none}"

    # Filter ClinVar to P/LP, >=2-star (and optionally rename contigs).
    python3 ${projectDir}/bin/filter_clinvar.py \\
        "${clinvar_vcf}" \\
        --out clinvar.filtered.vcf.gz \\
        \${RENAME_FLAG}

    # Cross-reference against sample. If CRAM is provided, also emit verification
    # shell scripts that reproduce the manual three-check workflow (VCF record,
    # pileup, MAPQ distribution).
    CRAM_ARGS=""
    if [[ -s "${sample_cram}" && -s "${reference_fasta}" ]]; then
        CRAM_ARGS="--cram ${sample_cram} --reference-fasta ${reference_fasta}"
    fi
    python3 ${projectDir}/bin/annotate_findings.py \\
        --sample "${sample_vcf}" \\
        --clinvar clinvar.filtered.vcf.gz \\
        --acmg-sf "${acmg_sf_yaml}" \\
        --out-acmg-sf 02_acmg_sf_findings.tsv \\
        --out-other  03_clinvar_other_findings.tsv \\
        \${CRAM_ARGS}
    """
}
