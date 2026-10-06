/*
 * Extract specific headers from Bruker dot d folders
 **/
process BRUKERMETRICSEXTRACTION {
    tag "${meta.id}"
    label 'process_low'

    stageInMode 'copy'  // needed due to pyopenms

    container "ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-604eac1"

    input:
    tuple val(meta), path(dotd_bruker_folder)

    output:
    tuple val(meta), path("*.hdf5"), emit: hdf5
    tuple val("${task.process}"), val('macproqc_helpers'), eval('macproqc-helpers --version | sed "s;^macproqc-helpers ;;"'), topic: versions, emit: versions_macproqc_helpers

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.bruker_metrics"
    // alphatims.bruker.TimsTOF (used in macproqc_helpers) dispatches on the folder name's
    // suffix, requiring it to end in ".d". the staged input folder may not (e.g.
    // when staged under a pipeline sample ID), so alias it to a ".d"-suffixed name.
    def orig_name = dotd_bruker_folder.name.toString()
    def dotd_name = orig_name.endsWith('.d') ? orig_name : "${orig_name}.d"
    def link_cmd = dotd_name == orig_name ? '' : "ln -s ${dotd_bruker_folder} ${dotd_name}"
    // list values may be given as a list or as a comma-separated string
    def metric_args = [
        meta.bruker_headers_to_parse ? [ meta.bruker_headers_to_parse ].flatten().collectMany { elem -> elem.toString().split(',') as List }*.trim().findAll { elem -> elem }.collect { elem -> "-headers_to_parse '${elem}'" }.join(' ') : '',
        meta.bruker_frame_headers_to_parse ? [ meta.bruker_frame_headers_to_parse ].flatten().collectMany { elem -> elem.toString().split(',') as List }*.trim().findAll { elem -> elem }.collect { elem -> "-frame_headers_to_parse '${elem}'" }.join(' ') : '',
        meta.bruker_calibrants ? [ meta.bruker_calibrants ].flatten().collectMany { elem -> elem.toString().split(',') as List }*.trim().findAll { elem -> elem }.collect { elem -> "-calibrants_to_retrieve '${elem}'" }.join(' ') : '',
        meta.bruker_calibrants_mz_tolerance ? "-calibrants_mz_tolerance ${meta.bruker_calibrants_mz_tolerance}" : '',
        meta.bruker_calibrants_mobility_tolerance ? "-calibrants_mobility_tolerance ${meta.bruker_calibrants_mobility_tolerance}" : '',
    ].join(' ')
    """
    ${link_cmd}
    python -m macproqc_helpers collect-metrics-from-bruker \\
        ${metric_args} \\
        ${args} \\
        -d_folder ${dotd_name} \\
        -out_hdf5 ${prefix}.hdf5
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.bruker_metrics"
    """
    echo ${args}

    touch ${prefix}.hdf5
    """
}
