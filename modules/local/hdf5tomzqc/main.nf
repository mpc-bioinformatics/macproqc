process HDF5TOMZQC {
    tag "$meta.id"
    label 'process_single'

    container 'ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-604eac1'

    input:
    tuple val(meta), path(hdf5)

    output:
    tuple val(meta), path("${prefix}.mzqc"), emit: mzqc
    tuple val("${task.process}"), val('macproqc_helpers'), eval('macproqc-helpers --version | sed "s;^macproqc-helpers ;;"'), topic: versions, emit: versions_macproqc_helpers

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python -m macproqc_helpers hdf5-to-mzqc \\
        -hdf5 ${hdf5} \\
        -mzqc_out ${prefix}.mzqc \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # create an empty mzqc file
    touch ${prefix}.mzqc
    """
}
