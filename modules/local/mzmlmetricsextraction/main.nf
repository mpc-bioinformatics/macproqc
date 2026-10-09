/*
 * Extract specific metrics from mzML files
 **/
process MZMLMETRICSEXTRACTION {
    tag "${meta.id}"
    label 'process_low'

    container "ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-a9c31aa"

    input:
    tuple val(meta), path(mzml_file)

    output:
    tuple val(meta), path("*.hdf5"), emit: hdf5
    tuple val("${task.process}"), val('macproqc_helpers'), eval('macproqc-helpers --version | sed "s;^macproqc-helpers ;;"'), topic: versions, emit: versions_macproqc_helpers

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // the first two are required by the CLI, hence the fallbacks
    def metric_args = [
        "-report_up_to_charge ${meta.report_up_to_charge ?: 5}",
        meta.base_peak_tic_up_to ? "-base_peak_tic_up_to ${meta.base_peak_tic_up_to}" : '',
        meta.ms1_map_rt_bins ? "-ms1_map_rt_bins ${meta.ms1_map_rt_bins}" : '',
        meta.ms1_map_mz_bins ? "-ms1_map_mz_bins ${meta.ms1_map_mz_bins}" : '',
    ].join(' ')

    """
    python -m macproqc_helpers collect-metrics-from-mzml \\
        ${metric_args} \\
        ${args} \\
        -mzml ${mzml_file} \\
        -out_hdf5 ${prefix}.mzml_metrics.hdf5
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.mzml_metrics.hdf5

    """
}
