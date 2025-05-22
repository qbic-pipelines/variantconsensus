/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Used Modules
include { BCFTOOLS_VIEW as FILTER_SNPS   } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as FILTER_INDELS } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_ISEC as ISEC_SNPS     } from '../modules/nf-core/bcftools/isec/main'
include { TABIX_BGZIPTABIX as TABIX_SNPS } from '../modules/nf-core/tabix/bgziptabix/main'
include { BCFTOOLS_VIEW as PASS_SNPS     } from '../modules/nf-core/bcftools/view/main'

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

    ch_both = ch_samplesheet.filter { meta, _files ->
            meta.varianttype != 'snps' && meta.varianttype != 'indels'
        }

    ch_snps = ch_samplesheet.filter { meta, _files ->
            meta.varianttype == 'snps'
        }

    // TODO: Uncomment indel
    // ch_indels = ch_samplesheet.filter { meta, _files ->
    //         meta.varianttype == 'indels'
    //     }

    // FILTER_SNPS to divide into SNPs and INDELs
    FILTER_SNPS(
        ch_both.map { meta, vcfs -> [meta, vcfs[0], vcfs[1]] },
        [],
        [],
        [],
    )

    ch_divided_snps = FILTER_SNPS.out.vcf.join(FILTER_SNPS.out.tbi)
            .map { meta, vcf, tbi -> [ meta.subMap(meta.keySet() - ['varianttype']) + [ 'varianttype': 'snps' ], [vcf, tbi] ] }

    ch_versions = ch_versions.mix(FILTER_SNPS.out.versions)

    // TODO: uncomment INDEL filter
    // // FILTER_INDELS to divide into SNPs and INDELs
    // FILTER_INDELS(
    //     ch_both.map { meta, vcfs -> [meta, vcfs[0], vcfs[1]] },
    //     [],
    //     [],
    //     [],
    // )

    // ch_all_indels = ch_indels.mix(
    //     FILTER_INDELS.out.vcf.join(FILTER_INDELS.out.tbi)
    //         .map { meta, vcf, tbi -> [ meta.subMap(meta.keySet() - ['varianttype']) + [ 'varianttype': 'indels' ], [vcf, tbi] ] }
    //     )

    // ch_versions = ch_versions.mix(FILTER_INDELS.out.versions)

    // Define the consensus thresholds you want to compare
    def consensusThresholds = [1, 2, 3] // Adjust these values as needed

    // Group SNPs for ISEC
    ch_snps_grouped = Channel.empty()
        .mix(ch_snps, ch_divided_snps)
        .map { meta, files ->
            [ [meta.id, meta.sample, meta.varianttype], meta, files[0], files[1] ]}
        .groupTuple(by: [0])
        .map { _id, metas, vcfs, tbis ->
            def meta = metas[0].subMap(metas[0].keySet() - ['caller'])
            meta = meta + [ 'numFiles': vcfs.flatten().size() ]
            [meta, vcfs.flatten(), tbis.flatten()]
        }

    // Create multiple channels with different consensusFiles values
    ch_consensus_channels = consensusThresholds.collect { threshold ->
        ch_snps_grouped.map { meta, vcfs, tbis ->
            def newMeta = meta + [ 'consensusFiles': meta.numFiles - threshold ]
            [newMeta, vcfs, tbis]
        }
    }

    // Combine all channels into a single channel
    ch_all_consensus = Channel.empty()
    ch_consensus_channels.each { channel ->
        ch_all_consensus = ch_all_consensus.mix(channel)
    }

    ISEC_SNPS( ch_all_consensus )

    ISEC_SNPS.out.results
        .map { meta, dir ->
            def new_filename = "${meta.patient}.${meta.id}.${meta.varianttype}.N-${meta.consensusFiles}.consensus.vcf"
            def copied_file = file("${dir}/${new_filename}")

            if (workflow.stubRun) {
                // For stub runs, just create an empty file
                copied_file.text = ''
            } else {
                // For actual runs, perform the copy operation
                def files = dir.listFiles()
                def original_file = files.find { it.name == '0000.vcf' }
                if (original_file) {
                    original_file.copyTo(copied_file)
                } else {
                    log.warn "File '0000.vcf' not found in directory ${dir}. Creating an empty file."
                    copied_file.text = ''
                }
            }

            return [meta, copied_file]
        }
        .set { ch_intersect_all_snps }

    ch_versions = ch_versions.mix(ISEC_SNPS.out.versions)

    // TODO: BCFTOOLS ISEC for INDEL consensus


    // Zip and index the SNP consensus VCFs
    TABIX_SNPS( ch_intersect_all_snps )

    ch_versions = ch_versions.mix(TABIX_SNPS.out.versions)

    // TODO: TABIX for INDEL consensus

    // Filter the SNPs for PASS variants
    PASS_SNPS( TABIX_SNPS.out.gz_tbi, [], [], [] )

    ch_versions = ch_versions.mix(PASS_SNPS.out.versions)

    // TODO: PASS for INDEL consensus


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
