process CLINVAR_FETCH {
    tag 'NCBI ClinVar (GRCh38)'

    // Persistent cache so reruns don't re-download the ~70 MB VCF.
    storeDir "${params.clinvar_cache_dir}"

    conda 'bioconda::tabix=1.11'

    output:
    tuple path('clinvar.vcf.gz'), path('clinvar.vcf.gz.tbi'), path('clinvar.vcf.gz.md5'), path('clinvar_release_date.txt')

    script:
    """
    set -euo pipefail
    BASE_URL="${params.clinvar_url}"
    wget -q --no-clobber "\${BASE_URL}" -O clinvar.vcf.gz
    wget -q --no-clobber "\${BASE_URL}.md5" -O clinvar.vcf.gz.md5
    wget -q --no-clobber "\${BASE_URL}.tbi" -O clinvar.vcf.gz.tbi

    # Verify integrity
    md5sum -c clinvar.vcf.gz.md5

    # Capture the release date from the VCF header for provenance
    zcat clinvar.vcf.gz | head -100 | grep -E "^##(fileDate|source)" > clinvar_release_date.txt
    """
}
