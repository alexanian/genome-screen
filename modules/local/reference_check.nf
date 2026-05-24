process REFERENCE_CHECK {
    tag "${vcf.baseName}"

    conda 'bioconda::pyyaml=6.0'

    input:
    path vcf
    path signatures
    val  expected_build

    output:
    path '00_reference_check.txt', emit: report
    path '00_reference_check.json', emit: json

    script:
    """
    python3 ${projectDir}/bin/check_reference.py \\
        "${vcf}" \\
        --signatures "${signatures}" \\
        --expected "${expected_build}" \\
        > 00_reference_check.txt

    python3 ${projectDir}/bin/check_reference.py \\
        "${vcf}" \\
        --signatures "${signatures}" \\
        --expected "${expected_build}" \\
        --json \\
        > 00_reference_check.json
    """
}
