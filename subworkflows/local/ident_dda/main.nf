include { OPENMS_DECOYDATABASE } from '../../../modules/nf-core/openms/decoydatabase/main'
include { COMETCONFIG; COMETCONFIG as COMETCONFIG_LABELLED } from '../../../modules/local/cometconfig/main'
include { COMET } from '../../../modules/nf-core/comet/main'
include { PIA_CONFIGURE } from '../../../modules/local/pia/configure/main'
include { PIA_COMPILEXML; PIA_COMPILEXML as PIA_COMPILEXML_PREFILTERED } from '../../../modules/local/pia/compilexml/main'
include { PIA_ANALYSIS; PIA_ANALYSIS as PIA_ANALYSIS_PREFILTER; PIA_ANALYSIS as PIA_ANALYSIS_LABELLED } from '../../../modules/local/pia/analysis/main'
include { IDMETRICSEXTRACTION } from '../../../modules/local/idmetricsextraction/main'

workflow IDENT_DDA {
    take:
    fasta                   // val: fasta path
    skip_decoy_generation   // val: true to skip the decoy generation
    comet_config_template   // val: path to comet config template (or null to use default)
    search_labels           // val: true to search for label modifications, false to skip (which label is defined in modules.config)
    ch_mzml                 // channel: [ val(meta), mzmls ]
    pia_fdr_threshold       // val: threshold for PIA FDR filtering
    pia_prefilter_threshold // val: threshold for PIA prefiltering, if 0 or null/empty, no prefiltering is done

    main:

    // for now, only one FASTA is used for all runs
    ch_fasta = channel.fromPath(fasta, checkIfExists: true)
        .map { f -> [[id: "global_fasta"], f] }

    if (!skip_decoy_generation) {
        OPENMS_DECOYDATABASE(
            ch_fasta
        )
        ch_fasta = OPENMS_DECOYDATABASE.out.decoy_fasta
    }

    // create fiel channel for comet config template, either user-specified or default
    ch_comet_config_template = comet_config_template
        ? channel.fromPath(comet_config_template, checkIfExists: true)
        : channel.fromPath("${projectDir}/assets/default_configs/comet.params", checkIfExists: true)

    // add an id to the comet config template channel
    ch_comet_config_template = ch_comet_config_template.map { params ->
        def params_id = comet_config_template
            ? 'user-template'
            : 'default'
        [ [id: params_id], params]
    }

    // create the config file (for unlabelled search)
    COMETCONFIG(
        ch_comet_config_template
    )

    // combine the channels for unlabelled search
    ch_comet_input = ch_mzml.combine(ch_fasta).combine(COMETCONFIG.out.params)
        .map { mzml_meta, mzml, fasta_meta, fastafile, params_meta, comet_params ->
            def meta = mzml_meta + [mzml_id: mzml_meta.id, label_search: false, fasta_id: fasta_meta.id, params_id: params_meta.id]
            meta.id = "${meta.mzml_id}-unlabelled"
            [meta, mzml, fastafile, comet_params]
        }

    if (search_labels == true) {
        // create the config file (for labelled search)
        COMETCONFIG_LABELLED(
            ch_comet_config_template
        )

        // create channels for labelled search
        ch_comet_input_labelled = ch_mzml.combine(ch_fasta).combine(COMETCONFIG_LABELLED.out.params)
            .map { mzml_meta, mzml, fasta_meta, fastafile, params_meta, comet_params ->
                def meta = mzml_meta + [mzml_id: mzml_meta.id, label_search: true, fasta_id: fasta_meta.id, params_id: params_meta.id]
                meta.id = "${meta.mzml_id}-labelled"
                [meta, mzml, fastafile, comet_params]
            }

        // mix all params for comet call
        ch_comet_input = ch_comet_input.mix(ch_comet_input_labelled)
    }

    // perform identification
    COMET(
        ch_comet_input
    )

    ///////////////////////////////////////////////////////
    // PIA analysis

    // first round of PIA compilation (pre-filtereing or branching off labelled is performed later)
    PIA_COMPILEXML(
        COMET.out.mzid
    )

    // Create PIA configuration files:
    // - the main (final) analysis
    // - an optional pre-filter pass (for the unlabelled search)
    // - if labelled searches are performed, a PSM-only pass for those
    ch_pia_configuration_in = channel.of(
        [
            id: "main-analysis",
            psm_export: true,
            peptide_export: true,
            protein_export: true,
            fdr_filter: true,
            remove_decoys: true,
            fdr_threshold: pia_fdr_threshold,
            psm_export_format: 'mzTab'
        ]
    )

    if (pia_prefilter_threshold != null && pia_prefilter_threshold > 0) {
        ch_pia_configuration_in = ch_pia_configuration_in.mix(
            channel.of([
                id: "prefilter",
                psm_export: true,
                peptide_export: false,
                protein_export: false,
                fdr_filter: true,
                remove_decoys: false,
                fdr_threshold: pia_prefilter_threshold,
                psm_export_format: 'mzid'           // pre-filtering uses mzid, as it is more complete than mzTab
            ])
        )
    }

    if (search_labels == true) {
        ch_pia_configuration_in = ch_pia_configuration_in.mix(
            channel.of([
                id: "labelled-analysis",
                psm_export: true,
                peptide_export: false,
                protein_export: false,
                fdr_filter: false,
                remove_decoys: true,
                fdr_threshold: pia_fdr_threshold,
                psm_export_format: 'mzTab'
            ])
        )
    }

    PIA_CONFIGURE(ch_pia_configuration_in)


    // split the compiled PIA xml files for labelled and unlabelled searches
    ch_pia_xml_in = PIA_COMPILEXML.out.pia_xml.branch { meta, _xml ->
        unlabelled: meta.label_search == false
        labelled: meta.label_search == true
    }

    // for unlabelled: optionally pre-filter on PSM level, recompile, then run the final analysis
    ch_unlabelled_xml_for_main = ch_pia_xml_in.unlabelled

    if (pia_prefilter_threshold != null && pia_prefilter_threshold > 0) {
        // perform the pre-filtering and update the channel for teh main PIA analysis
        ch_prefilter_json = PIA_CONFIGURE.out.pia_json
            .filter { meta, _json -> meta.id == "prefilter" }
            .map { _meta, json -> json }

        PIA_ANALYSIS_PREFILTER(
            ch_pia_xml_in.unlabelled.combine(ch_prefilter_json)
        )

        PIA_COMPILEXML_PREFILTERED(
            PIA_ANALYSIS_PREFILTER.out.psms
        )

        ch_unlabelled_xml_for_main = PIA_COMPILEXML_PREFILTERED.out.pia_xml
    }

    // pick out the main analysis configuration JSON
    ch_main_json = PIA_CONFIGURE.out.pia_json
        .filter { meta, _json -> meta.id == "main-analysis" }
        .map { _meta, json -> json }

    PIA_ANALYSIS(
        ch_unlabelled_xml_for_main.combine(ch_main_json)
    )

    // for labelled: perform the PIA analysis only on PSM level, no pre-filtering
    ch_labelled_psms = channel.empty()

    if (search_labels == true) {
        ch_labelled_json = PIA_CONFIGURE.out.pia_json
            .filter { meta, _json -> meta.id == "labelled-analysis" }
            .map { _meta, json -> json }

        PIA_ANALYSIS_LABELLED(
            ch_pia_xml_in.labelled.combine(ch_labelled_json)
        )

        ch_labelled_psms = PIA_ANALYSIS_LABELLED.out.psms
    }

    // extract the identification metrics from the PIA results (unlabelled only, labelled will be used for XIC extraction in the respective module)
    ch_idmetrics_extraction_in = PIA_ANALYSIS.out.psms
        .join(PIA_ANALYSIS.out.peptides)
        .join(PIA_ANALYSIS.out.proteins)

    IDMETRICSEXTRACTION(
        ch_idmetrics_extraction_in
    )

    emit:
    mzid         = COMET.out.mzid
    pia_psms     = PIA_ANALYSIS.out.psms.mix(ch_labelled_psms)
    pia_peptides = PIA_ANALYSIS.out.peptides
    pia_proteins = PIA_ANALYSIS.out.proteins
    idmetrics    = IDMETRICSEXTRACTION.out.hdf5
}
