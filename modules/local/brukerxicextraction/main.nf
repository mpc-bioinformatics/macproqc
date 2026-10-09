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
    // alphatims.bruker.TimsTOF (used in macproqc_helpers) dispatches on the folder name's
    // suffix, requiring it to end in ".d". the staged input folder may not (e.g.
    // when staged under a pipeline sample ID), so alias it to a ".d"-suffixed name.
    def orig_name = d_folder.name.toString()
    def dotd_name = orig_name.endsWith('.d') ? orig_name : "${orig_name}.d"
    def link_cmd = dotd_name == orig_name ? '' : "ln -s ${d_folder} ${dotd_name}"

    """
    ${link_cmd}
    python -m macproqc_helpers extract-xic-bruker \\
        ${args} \\
        -d_folder ${dotd_name} \\
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
