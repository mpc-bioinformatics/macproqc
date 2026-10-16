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
include { COLLECT_METRICS } from '../subworkflows/local/collect_metrics'

include { HDF5TOMZQC } from '../modules/local/hdf5tomzqc'
include { QCVISUALIZATION } from '../modules/local/qcvisualization'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MACPROQC {

    take:
    ch_samplesheet              // channel: samplesheet read in from --input
    outdir
    spike_ins_table             // val: path to the spike-ins reference table (csv), or null to skip any spike in metric extraction
    fasta
    skip_decoy_generation
    comet_config_template       // val: path to comet config template, or null to use created default config
    label_modifications
    pia_fdr_threshold           // number: pia fdr threshold (default 0.01)
    pia_prefilter_threshold     // number: pia prefilter threshold (the FDR pre-filter or 0 / empty / null to skip prefiltering)
    output_table_type           // val: output table type for the QC visualization (default: csv)
    figure_format               // val: figure format for the QC visualization (default: plotly)

    main:

    def ch_versions = channel.empty()   // this channel needs to be filled, if the module/subworkflow return any versions NOT in the topic channel, but in the versions.yml file

    // prepare the spectra files (strip the fasta_file before passing to PREPARE_SPECTRA)
    PREPARE_SPECTRA(
        ch_samplesheet.map { meta, spectrum_file -> [meta, spectrum_file] }
    )

    // enable "label search" if the user has specified any label modifications (see also modules.config)
    def search_label_modifications = label_modifications?.trim() ? true : false

    // perform DDA identification and extract the respective metrics
    IDENT_DDA(
        fasta,
        skip_decoy_generation,
        comet_config_template,
        search_label_modifications,
        PREPARE_SPECTRA.out.mzmls,
        pia_fdr_threshold,
        pia_prefilter_threshold
    )

    // collect all metrics (from mzML, raw, and PIA results) and optionally also the spike-in XIC metrics
    COLLECT_METRICS(
        PREPARE_SPECTRA.out.mzmls,
        PREPARE_SPECTRA.out.uncompressed,
        IDENT_DDA.out.pia_psms
            .join(IDENT_DDA.out.pia_peptides)
            .join(IDENT_DDA.out.pia_proteins),
        IDENT_DDA.out.pia_labelled_psms,
        spike_ins_table
    )

    // export the metrics to mzQC format
    HDF5TOMZQC(
        COLLECT_METRICS.out.hdf5
    )

    // perform the visualization once over all runs (collect every per-run HDF5 into a single list, sorted by file name for a stable order)
    QCVISUALIZATION (
        COLLECT_METRICS.out.hdf5
            .map { _meta, hdf5 -> hdf5 }
            .collect(sort: { a, b -> a.name <=> b.name })
            .map { hdf5_files -> [ [ id: 'qc_visualization' ], hdf5_files ] },
        output_table_type,
        figure_format,
        spike_ins_table ? file(spike_ins_table, checkIfExists: true) : []   // [] = no spike-ins table
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
