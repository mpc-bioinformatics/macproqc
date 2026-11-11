/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_macproqc_pipeline'
include { warnMissingRuns        } from '../subworkflows/local/utils_nfcore_macproqc_pipeline'

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
    ch_samplesheet              // channel: samplesheet read in from --input (empty, if the metrics extraction is skipped)
    ch_input_hdf5               // channel: [ val(meta), path(hdf5) ] from --input_hdf5_dir (empty, if the metrics extraction is performed)
    skip_metrics_extraction     // boolean: skip the metrics extraction (spectra preparation, identification, metric collection) and use ch_input_hdf5
    skip_identification         // boolean: skip the identification (Comet + PIA) during the metrics extraction
    skip_mzqc_export            // boolean: skip the mzQC export
    skip_visualization          // boolean: skip the visualization
    ignore_failed_runs          // boolean: failed tasks reading the raw files are ignored, warn about the affected runs
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

    // the per-run metric HDF5 files: given via --input_hdf5_dir, or created by the metrics extraction
    def ch_hdf5 = ch_input_hdf5

    if (!skip_metrics_extraction) {
        // prepare the spectra files (strip the fasta_file before passing to PREPARE_SPECTRA)
        PREPARE_SPECTRA(
            ch_samplesheet.map { meta, spectrum_file -> [meta, spectrum_file] }
        )

        // without identification, no ID metrics are collected and the spike-in metrics use only the RT windows of the table
        def ch_pia_results = channel.empty()
        def ch_labelled_psm_mztab = channel.empty()

        if (!skip_identification) {
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

            ch_pia_results = IDENT_DDA.out.pia_psms
                .join(IDENT_DDA.out.pia_peptides)
                .join(IDENT_DDA.out.pia_proteins)
            ch_labelled_psm_mztab = IDENT_DDA.out.pia_labelled_psms
        }

        // collect all metrics (from mzML, raw, and PIA results) and optionally also the spike-in XIC metrics
        COLLECT_METRICS(
            PREPARE_SPECTRA.out.mzmls,
            PREPARE_SPECTRA.out.uncompressed,
            ch_pia_results,
            ch_labelled_psm_mztab,
            spike_ins_table
        )

        ch_hdf5 = COLLECT_METRICS.out.hdf5

        if (ignore_failed_runs) {
            // warn about runs with failed (and ignored) decompression / mzML conversion, or without any results
            warnMissingRuns(
                ch_samplesheet,
                PREPARE_SPECTRA.out.uncompressed,
                { meta, file -> "Run '${meta.id}' (${file}) could not be decompressed (or has an unsupported format)." }
            )
            warnMissingRuns(
                PREPARE_SPECTRA.out.uncompressed,
                PREPARE_SPECTRA.out.mzmls,
                { meta, file -> "Run '${meta.id}' (${file.name}) could not be converted to mzML: mzML and identification metrics are missing." }
            )
            // the vendor specific extractions only run on (uncompressed) Thermo/Bruker files
            def ch_vendor_raw = PREPARE_SPECTRA.out.uncompressed.filter { meta, _file -> meta.vendor in ['thermo', 'bruker'] }
            warnMissingRuns(
                ch_vendor_raw,
                COLLECT_METRICS.out.vendor_metrics,
                { meta, file -> "Run '${meta.id}' (${file.name}): the instrument header metrics extraction failed, these metrics are missing." }
            )
            if (spike_ins_table) {
                warnMissingRuns(
                    ch_vendor_raw,
                    COLLECT_METRICS.out.xic,
                    { meta, file -> "Run '${meta.id}' (${file.name}): the spike-in XIC extraction failed, spike-in metrics are missing." }
                )
            }
            warnMissingRuns(
                ch_samplesheet,
                ch_hdf5,
                { meta, file -> "!!! Run '${meta.id}' (${file}) could not be processed and is MISSING IN ALL RESULTS !!!" }
            )
        }
    }

    if (!skip_mzqc_export) {
        // export the metrics to mzQC format
        HDF5TOMZQC(
            ch_hdf5
        )
    }

    if (!skip_visualization) {
        // perform the visualization once over all runs (collect every per-run HDF5 into a single list, sorted by file name for a stable order)
        QCVISUALIZATION (
            ch_hdf5
                .map { _meta, hdf5 -> hdf5 }
                .collect(sort: { a, b -> a.name <=> b.name })
                .map { hdf5_files -> [ [ id: 'qc_visualization' ], hdf5_files ] },
            output_table_type,
            figure_format,
            spike_ins_table ? file(spike_ins_table, checkIfExists: true) : []   // [] = no spike-ins table
        )
    }

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
