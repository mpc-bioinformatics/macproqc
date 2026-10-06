/*
 * Extract specific metrics from mzML files
 **/
process MZMLMETRICSEXTRACTION {
    tag "${meta.id}"
    label 'process_low'

    container "ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-604eac1"

    input:
    tuple val(meta), path(mzml_file)

    output:
    tuple val(meta), path("*.hdf5"), emit: hdf5
    tuple val("${task.process}"), val('macproqc_helpers'), eval('macproqc-helpers --version | sed "s;^macproqc-helpers ;;"'), topic: versions, emit: versions_macproqc_helpers

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.mzml_metrics"
    // the first two are required by the CLI, hence the fallbacks
    def metric_args = [
        "-base_peak_tic_up_to ${meta.base_peak_tic_up_to ?: 9999}",
        "-report_up_to_charge ${meta.report_up_to_charge ?: 5}",
        meta.ms1_map_rt_bins ? "-ms1_map_rt_bins ${meta.ms1_map_rt_bins}" : '',
        meta.ms1_map_mz_bins ? "-ms1_map_mz_bins ${meta.ms1_map_mz_bins}" : '',
    ].join(' ')

    """
    python -m macproqc_helpers collect-metrics-from-mzml \\
        ${metric_args} \\
        ${args} \\
        -mzml ${mzml_file} \\
        -out_hdf5 ${prefix}.hdf5
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.mzml_metrics"
    """
    echo ${args}

    touch ${prefix}.hdf5

    """
}
