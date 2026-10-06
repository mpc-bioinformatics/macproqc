/*
 * Extract specific headers from Thermofisher RAW files
 **/
process THERMOMETRICSEXTRACTION {
    tag "${meta.id}"
    label 'process_low'

    stageInMode 'copy'  // needed due to mono

    container "ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-604eac1"

    input:
    tuple val(meta), file(raw_thermo_file)

    output:
    tuple val(meta), path("*.hdf5"), emit: hdf5
    tuple val("${task.process}"), val('macproqc_helpers'), eval('macproqc-helpers --version | sed "s;^macproqc-helpers ;;"'), topic: versions, emit: versions_macproqc_helpers

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.thermo_metrics"
    // list values may be given as a list or as a comma-separated string
    def metric_args = [
        meta.thermo_extra_headers_to_parse ? [ meta.thermo_extra_headers_to_parse ].flatten().collectMany { elem -> elem.toString().split(',') as List }*.trim().findAll { elem -> elem }.collect { elem -> "-extra_headers_to_parse '${elem}'" }.join(' ') : '',
        meta.thermo_tune_headers_to_parse ? [ meta.thermo_tune_headers_to_parse ].flatten().collectMany { elem -> elem.toString().split(',') as List }*.trim().findAll { elem -> elem }.collect { elem -> "-tune_headers_to_parse '${elem}'" }.join(' ') : '',
        meta.thermo_log_headers_to_parse ? [ meta.thermo_log_headers_to_parse ].flatten().collectMany { elem -> elem.toString().split(',') as List }*.trim().findAll { elem -> elem }.collect { elem -> "-log_headers_to_parse '${elem}'" }.join(' ') : '',
    ].join(' ')
    """
    python -m macproqc_helpers collect-metrics-from-thermo \\
        ${metric_args} \\
        ${args} \\
        -raw ${raw_thermo_file} \\
        -out_hdf5 ${prefix}.hdf5
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.thermo_metrics"
    """
    echo ${args}

    touch ${prefix}.hdf5
    """
}
