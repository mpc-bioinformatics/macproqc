/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_macproqc_pipeline'

include { PREPARE_SPECTRA } from '../subworkflows/local/prepare_spectra'
include { IDENT_DDA } from '../subworkflows/local/ident_dda'

include { COMBINEHDF5 } from '../modules/local/combinehdf5/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MACPROQC {

    take:
    ch_samplesheet              // channel: samplesheet read in from --input
    outdir
    fasta
    skip_decoy_generation
    comet_config_template       // string: path to comet config template, or null to use created default config
    label_modifications         //
    pia_fdr_threshold           // number: pia fdr threshold (default 0.01)
    pia_prefilter_threshold     // number: pia prefilter threshold (the FDR pre-filter or 0 / empty / null to skip prefiltering)

    main:

    def ch_versions = channel.empty()   // this channel needs to be filled, if the module/subworkflow return any versions NOT in the topic channel, but in the versions.yml file

    // prepare the spectra files (strip the fasta_file before passing to PREPARE_SPECTRA)
    PREPARE_SPECTRA(
        ch_samplesheet.map { meta, spectrum_file -> [meta, spectrum_file] }
    )

    // enable "label search" if the user has specified any label modifications (see also modules.config)
    def search_label_modifications = label_modifications?.trim() ? true : false

    // perform DDA identification nad extract the respective metrics
    IDENT_DDA(
        fasta,
        skip_decoy_generation,
        comet_config_template,
        search_label_modifications,
        PREPARE_SPECTRA.out.mzmls,
        pia_fdr_threshold,
        pia_prefilter_threshold
    )

    // combine the  metrics from all analyses runs into one HDF5 file
    COMBINEHDF5(
        IDENT_DDA.out.idmetrics
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'macproqc_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )
    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
