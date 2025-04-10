/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Used Modules
include { BCFTOOLS_VIEW as FILTER_SNPS   } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as FILTER_INDELS } from '../modules/nf-core/bcftools/view/main'

// Template Modules
include { MULTIQC                        } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap               } from 'plugin/nf-schema'
include { paramsSummaryMultiqc           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML         } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText         } from '../subworkflows/local/utils_nfcore_variantconsensus_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANTCONSENSUS {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Check if the input contains both SNPs and INDELs
    ch_snps = ch_samplesheet.filter { meta, vcfs ->
        meta.containsKey('varianttype') && meta.varianttype == 'snps'
    }.map { meta, vcfs -> [meta, vcfs[0], vcfs[1]] }

    ch_indels = ch_samplesheet.filter { meta, vcfs ->
        meta.containsKey('varianttype') && meta.varianttype == 'indels'
    }.map { meta, vcfs -> [meta, vcfs[0], vcfs[1]] }

    ch_unknown = ch_samplesheet.filter { meta, vcfs ->
        meta.containsKey('varianttype') && meta.varianttype == 'both'
    }

    // FILTER_SNPS to divide into SNPs and INDELs
    FILTER_SNPS(
        ch_unknown.map { meta, vcfs -> [meta, vcfs[0], vcfs[1]] },
        [],
        [],
        [],
    )

    ch_all_snps = ch_snps.mix(
        FILTER_SNPS.out.vcf.join(FILTER_SNPS.out.tbi)
            .map { meta, vcf, tbi -> [ meta.subMap(meta.keySet() - ['varianttype']) + [ 'varianttype': 'snps' ], vcf, tbi ] }
        )

    ch_versions = ch_versions.mix(FILTER_SNPS.out.versions)

    // FILTER_INDELS to divide into SNPs and INDELs
    FILTER_INDELS(
        ch_unknown.map { meta, vcfs -> [meta, vcfs[0], vcfs[1]] },
        [],
        [],
        [],
    )

    ch_all_indels = ch_indels.mix(
        FILTER_INDELS.out.vcf.join(FILTER_INDELS.out.tbi)
            .map { meta, vcf, tbi -> [ meta.subMap(meta.keySet() - ['varianttype']) + [ 'varianttype': 'indels' ], vcf, tbi ] }
        )

    ch_versions = ch_versions.mix(FILTER_INDELS.out.versions)

    ch_all_snps.dump(tag: 'snps-filtered')
    ch_all_indels.dump(tag: 'indels-filtered')

    // TODO: BCFTOOLS ISEC for SNP consensus


    // TODO: BCFTOOLS ISEC for INDEL consensus


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'variantconsensus_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
