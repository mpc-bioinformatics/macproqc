/*
 * Collect spike-in/reference peptide metrics (maximum XIC intensity, retention time at
 * that maximum, number of identifications) from extracted XICs and identifications, and
 * store them in HDF5 format
 **/
process SPIKEINMETRICSEXTRACTION {
    tag "${meta.id}"
    label 'process_single'

    container "ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-604eac1"

    input:
    tuple val(meta), path(xic_json), path(identifications)
    path spike_ins_table

    output:
    tuple val(meta), path("*.hdf5"), emit: hdf5
    tuple val("${task.process}"), val('macproqc_helpers'), eval('macproqc-helpers --version | sed "s;^macproqc-helpers ;;"'), topic: versions, emit: versions_macproqc_helpers

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.spikein_metrics"
    """
    python -m macproqc_helpers collect-spikein-metrics \\
        ${args} \\
        -itrfp_json ${xic_json} \\
        -iidentifications ${identifications} \\
        -ispikeins ${spike_ins_table} \\
        -ohdf5 ${prefix}.hdf5
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.spikein_metrics"
    """
    echo ${args}

    touch ${prefix}.hdf5
    """
}
