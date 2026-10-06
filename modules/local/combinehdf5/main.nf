/*
 * Merge all per-metric HDF5 files produced for a single run into one combined HDF5 file
 **/
process COMBINEHDF5 {
    tag "${meta.id}"
    label 'process_single'

    container "ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-a9c31aa"

    input:
    tuple val(meta), path(hdf5_files)

    output:
    tuple val(meta), path("*.hdf5"), emit: hdf5
    tuple val("${task.process}"), val('macproqc_helpers'), eval('macproqc-helpers --version | sed "s;^macproqc-helpers ;;"'), topic: versions, emit: versions_macproqc_helpers

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def metric_args = [
        meta.hdf5_put_under_subdataset ? '-put_under_subdataset' : '',
        meta.hdf5_write_metadata ? '-write_metadata' : '',
    ].join(' ')
    """
    python -m macproqc_helpers combine-hdf5 \\
        ${metric_args} \\
        ${args} \\
        -hdf_out_name ${prefix}.hdf5 \\
        ${hdf5_files}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.hdf5
    """
}
