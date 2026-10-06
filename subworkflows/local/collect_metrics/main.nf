include { MZMLMETRICSEXTRACTION }    from '../../../modules/local/mzmlmetricsextraction/main'
include { THERMOMETRICSEXTRACTION }  from '../../../modules/local/thermometricsextraction/main'
include { BRUKERMETRICSEXTRACTION }  from '../../../modules/local/brukermetricsextraction/main'
include { XICEXTRACTIONCONFIG }      from '../../../modules/local/xicextractionconfig/main'
include { THERMOXICEXTRACTION }      from '../../../modules/local/thermoxicextraction/main'
include { BRUKERXICEXTRACTION }      from '../../../modules/local/brukerxicextraction/main'
include { SPIKEINMETRICSEXTRACTION } from '../../../modules/local/spikeinmetricsextraction/main'
include { COMBINEHDF5 }              from '../../../modules/local/combinehdf5/main'

// optional extraction settings can be set per sample in the samplesheet; samples without a value fall back to the pipeline param of the same name
def withMetricSettings(meta) {
    def keys = [
        'base_peak_tic_up_to', 'report_up_to_charge', 'ms1_map_rt_bins', 'ms1_map_mz_bins',
        'bruker_headers_to_parse', 'bruker_frame_headers_to_parse', 'bruker_calibrants',
        'bruker_calibrants_mz_tolerance', 'bruker_calibrants_mobility_tolerance',
        'thermo_extra_headers_to_parse', 'thermo_tune_headers_to_parse', 'thermo_log_headers_to_parse',
        'hdf5_put_under_subdataset', 'hdf5_write_metadata',
    ]
    def settings = keys.collectEntries { key -> [ (key): meta[key] != null ? meta[key] : params[key] ] }.findAll { _key, value -> value != null }
    return meta + settings
}

workflow COLLECT_METRICS {
    take:
    ch_mzml            // channel: [ val(meta), path(mzml) ]
    ch_raw             // channel: [ val(meta), path(raw|d_folder) ] with meta.vendor in ['thermo', 'bruker']
    ch_psm_mztab       // channel: [ val(meta), path(psm_mztab_file) ] — only consumed when collect_spike_ins is true
    spike_ins_table    // val: path to the spike-ins reference table (csv) — only consumed when collect_spike_ins is true
    collect_spike_ins  // val: true to also run the spike-in XIC metrics chain

    main:

    // add the optional extraction settings to the meta maps (samplesheet value, else pipeline param)
    ch_mzml = ch_mzml.map { meta, mzml -> [ withMetricSettings(meta), mzml ] }
    ch_raw  = ch_raw.map { meta, raw -> [ withMetricSettings(meta), raw ] }

    // mzML-based metrics (MS1/MS2 counts, TIC, RT, precursor charge, ...)
    MZMLMETRICSEXTRACTION(ch_mzml)

    // vendor-specific instrument header metrics (pump pressure, calibrants, ...)
    ch_raw
        .branch { meta, _raw ->
            thermo: meta.vendor == 'thermo'
            bruker: meta.vendor == 'bruker'
        }
        .set { ch_branched_raw }

    THERMOMETRICSEXTRACTION(ch_branched_raw.thermo)
    BRUKERMETRICSEXTRACTION(ch_branched_raw.bruker)

    ch_all_hdf5 = MZMLMETRICSEXTRACTION.out.hdf5
        .mix(THERMOMETRICSEXTRACTION.out.hdf5)
        .mix(BRUKERMETRICSEXTRACTION.out.hdf5)

    if (collect_spike_ins) {
        ch_spike_ins_table = channel.fromPath(spike_ins_table, checkIfExists: true)

        // create the XIC extraction config and the identifications matched to the spike-ins, based on the PSM results
        XICEXTRACTIONCONFIG(
            ch_psm_mztab,
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
