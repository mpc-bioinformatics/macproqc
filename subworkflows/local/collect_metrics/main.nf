include { IDMETRICSEXTRACTION }    from '../../../modules/local/idmetricsextraction/main'
include { MZMLMETRICSEXTRACTION }    from '../../../modules/local/mzmlmetricsextraction/main'
include { THERMOMETRICSEXTRACTION }  from '../../../modules/local/thermometricsextraction/main'
include { BRUKERMETRICSEXTRACTION }  from '../../../modules/local/brukermetricsextraction/main'
include { XICEXTRACTIONCONFIG }      from '../../../modules/local/xicextractionconfig/main'
include { THERMOXICEXTRACTION }      from '../../../modules/local/thermoxicextraction/main'
include { BRUKERXICEXTRACTION }      from '../../../modules/local/brukerxicextraction/main'
include { SPIKEINMETRICSEXTRACTION } from '../../../modules/local/spikeinmetricsextraction/main'
include { COMBINEHDF5 }              from '../../../modules/local/combinehdf5/main'

workflow COLLECT_METRICS {
    take:
    ch_mzml                     // channel: [ val(meta), path(mzml) ]
    ch_raw                      // channel: [ val(meta), path(raw|d_folder) ] with meta.vendor in ['thermo', 'bruker']
    ch_unlabelled_pia_results   // channel: [ val(meta), path(pia_psms), path(pia_peptides), path(pia_proteins) ] - the (unlabelled) ID results
    ch_labelled_psm_mztab       // channel: [ val(meta), path(psm_mztab_file) ] — only consumed when spike_ins_table is given
    spike_ins_table             // val: path to the spike-ins reference table (csv), or null to skip the spike-in XIC metrics chain

    main:

    // this channel will be filled with all per-metric HDF5 produced for a given run, to be merged into one HDF5 at the end
    ch_all_hdf5 = channel.empty()

    // Perform ID based metric extraction (if identification data is provided)
    if (ch_unlabelled_pia_results) {
        IDMETRICSEXTRACTION(ch_unlabelled_pia_results)
        ch_all_hdf5 = ch_all_hdf5.mix(IDMETRICSEXTRACTION.out.hdf5)
    }

    // mzML-based metrics (MS1/MS2 counts, TIC, RT, precursor charge, ...)
    // TODO: check the meta.params for the MZMLMETRICSEXTRACTION, nothing is forwarded right now! Add these as params!
    MZMLMETRICSEXTRACTION(ch_mzml)
    ch_all_hdf5 = ch_all_hdf5.mix(MZMLMETRICSEXTRACTION.out.hdf5)

    // vendor-specific instrument header metrics (pump pressure, calibrants, ...)
    ch_raw
        .branch { meta, _raw ->
            thermo: meta.vendor == 'thermo'
            bruker: meta.vendor == 'bruker'
        }
        .set { ch_branched_raw }

    // TODO: check the meta.params for the THERMOMETRICSEXTRACTION, nothing is forwarded right now! Add these as params!
    THERMOMETRICSEXTRACTION(ch_branched_raw.thermo)
    // TODO: check the meta.params for the BRUKERMETRICSEXTRACTION, nothing is forwarded right now! Add these as params!
    BRUKERMETRICSEXTRACTION(ch_branched_raw.bruker)
    ch_all_hdf5 = ch_all_hdf5
        .mix(THERMOMETRICSEXTRACTION.out.hdf5)
        .mix(BRUKERMETRICSEXTRACTION.out.hdf5)

    if (spike_ins_table) {
        ch_spike_ins_table = channel.value(file(spike_ins_table, checkIfExists: true))

        // create the XIC extraction config and the identifications matched to the spike-ins, based on the PSM results
        XICEXTRACTIONCONFIG(
            ch_labelled_psm_mztab,
            ch_spike_ins_table
        )

        // pair each raw file with its XIC extraction config (joined on meta.id, since the meta maps
        // themselves may carry different extra keys depending on their origin)
        ch_raw_and_config = ch_raw
            .map { meta, raw -> [ meta.id, meta, raw ] }
            .join(
                XICEXTRACTIONCONFIG.out.config.map { meta, xic_config, ident_csv -> [ meta.id, xic_config, ident_csv ] },
                by: 0
            )
            .map { _id, meta, raw, xic_config, ident_csv -> [ meta, raw, xic_config, ident_csv ] }

        // branch by vendor to run the correct XIC extraction (mzML-only runs have no raw/d_folder and are dropped here)
        ch_raw_and_config
            .branch { meta, _raw, _xic_config, _ident_csv ->
                thermo: meta.vendor == 'thermo'
                bruker: meta.vendor == 'bruker'
            }
            .set { ch_branched_raw_and_config }

        THERMOXICEXTRACTION(
            ch_branched_raw_and_config.thermo.map { meta, raw, xic_config, _ident_csv -> [ meta, raw, xic_config ] }
        )
        BRUKERXICEXTRACTION(
            ch_branched_raw_and_config.bruker.map { meta, raw, xic_config, _ident_csv -> [ meta, raw, xic_config ] }
        )

        ch_xics = THERMOXICEXTRACTION.out.xic.mix(BRUKERXICEXTRACTION.out.xic)

        // re-attach the identifications csv (joined on meta.id) for the metrics extraction
        ch_xics_and_identifications = ch_xics
            .map { meta, xic -> [ meta.id, meta, xic ] }
            .join(
                ch_raw_and_config.map { meta, _raw, _xic_config, ident_csv -> [ meta.id, ident_csv ] },
                by: 0
            )
            .map { _id, meta, xic, ident_csv -> [ meta, xic, ident_csv ] }

        SPIKEINMETRICSEXTRACTION(
            ch_xics_and_identifications,
            ch_spike_ins_table
        )

        ch_all_hdf5 = ch_all_hdf5.mix(SPIKEINMETRICSEXTRACTION.out.hdf5)
    }

    // group every per-metric HDF5 produced for a given run and merge them into one HDF5
    // (metas may carry different keys, so merge all of them instead of picking one)
    ch_grouped_hdf5 = ch_all_hdf5
        .map { meta, hdf5 -> [ meta.id, meta, hdf5 ] }
        .groupTuple(by: 0)
        .map { _id, metas, hdf5s -> [ metas.inject([:]) { merged, m -> merged + m }, hdf5s ] }

    COMBINEHDF5(ch_grouped_hdf5)

    emit:
    hdf5 = COMBINEHDF5.out.hdf5   // channel: [ val(meta), path(hdf5) ]
}
