#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Module imports — most are stubs while the pipeline is built incrementally
include { REFERENCE_CHECK   } from './modules/local/reference_check.nf'
// include { VCF_QC            } from './modules/local/vcf_qc.nf'
// include { CLINVAR_FETCH     } from './modules/local/clinvar_fetch.nf'
// include { CLINVAR_INTERSECT } from './modules/local/clinvar_intersect.nf'
// include { CURATED_LOOKUP    } from './modules/local/curated_lookup.nf'
// include { NUTRIGENOMICS     } from './modules/local/nutrigenomics.nf'
// include { PHARMCAT          } from './modules/local/pharmcat.nf'
// include { REPORT            } from './modules/local/report.nf'

workflow {
    // Required inputs
    if (!params.sample_vcf)      error 'Missing required param: --sample_vcf'
    if (!params.reference_fasta) error 'Missing required param: --reference_fasta'

    sample_vcf_ch       = Channel.fromPath(params.sample_vcf,      checkIfExists: true)
    reference_sigs_ch   = Channel.fromPath(params.reference_signatures, checkIfExists: true)

    REFERENCE_CHECK(sample_vcf_ch, reference_sigs_ch, params.expected_build)

    // Downstream steps will be added as their modules land.
    // VCF_QC(sample_vcf_ch)
    // CLINVAR_FETCH() | CLINVAR_INTERSECT(sample_vcf_ch, params.acmg_sf_yaml)
    // CURATED_LOOKUP(sample_vcf_ch, sample_cram_ch, params.curated_risk_yaml)
    // NUTRIGENOMICS(sample_vcf_ch, sample_cram_ch, params.nutrigenomic_yaml)
    // PHARMCAT(sample_vcf_ch)
    // REPORT(...)
}
