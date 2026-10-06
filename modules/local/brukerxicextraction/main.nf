/*
 * Extract XICs for spike-in/reference peptides from Bruker .d-folders, based on a
 * ThermoRawFileParser-compatible XIC extraction configuration (see XICEXTRACTIONCONFIG)
 **/
process BRUKERXICEXTRACTION {
    tag "${meta.id}"
    label 'process_medium'

    stageInMode 'copy'  // needed due to alphatims

    container "ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-a9c31aa"

    input:
    tuple val(meta), path(d_folder), path(xic_config)

    output:
    tuple val(meta), path("*.json"), emit: xic
    tuple val("${task.process}"), val('macproqc_helpers'), eval('macproqc-helpers --version | sed "s;^macproqc-helpers ;;"'), topic: versions, emit: versions_macproqc_helpers

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python -m macproqc_helpers extract-xic-bruker \\
        ${args} \\
        -d_folder ${d_folder} \\
        -in_json ${xic_config} \\
        -out_json ${prefix}.json
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    echo '{}' > ${prefix}.json
    """
}
