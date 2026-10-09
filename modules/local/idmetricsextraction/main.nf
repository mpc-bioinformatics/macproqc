process IDMETRICSEXTRACTION {
    tag "$meta.id"
    label 'process_single'

    container 'ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-604eac1'

    input:
    tuple val(meta), path(psms), path(peptides), path(proteins)

    output:
    tuple val(meta), path("${prefix}.idmetrics.hdf5"), emit: hdf5
    tuple val("${task.process}"), val('macproqc_helpers'), eval('macproqc-helpers --version | sed "s;^macproqc-helpers ;;"'), topic: versions, emit: versions_macproqc_helpers

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python -m macproqc_helpers collect-pia-metrics \\
        --pia_PSMs ${psms} \\
        --pia_peptides ${peptides} \\
        --pia_proteins ${proteins} \\
        --out_hdf5 ${prefix}.idmetrics.hdf5 \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.idmetrics.hdf5
    """
}
