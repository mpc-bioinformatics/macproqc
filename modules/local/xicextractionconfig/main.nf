/*
 * Create a ThermoRawFileParser (also alphatims-compatible) XIC extraction config
 * from a spike-ins CSV and (optionally) a PSM mzTab file to refine the RT windows
 **/
process XICEXTRACTIONCONFIG {
    tag "${meta.id}"
    label 'process_single'

    container "ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-55b549b"

    input:
    tuple val(meta), path(psm_mztab_files)   // psm_mztab_files: one or more PSM mzTab files, or [] to use only the RT windows of the spike-ins table
    path spike_ins

    output:
    tuple val(meta), path("*.xic_config.json"), path("*.ident.csv"), emit: config   // distinct suffix, so the XIC extraction outputs ("<id>.json") do not clash with this input
    tuple val("${task.process}"), val('macproqc_helpers'), eval('macproqc-helpers --version | sed "s;^macproqc-helpers ;;"'), topic: versions, emit: versions_macproqc_helpers

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def idents_arg = psm_mztab_files.collect { mztab -> "-mztab ${mztab}" }.join(' ')   // empty for []

    """
    python -m macproqc_helpers create_thermorawfileparser_xic_config \\
        ${args} \\
        -spikeins ${spike_ins} \\
        -xic_definition_json_out ${prefix}.xic_config.json \\
        -id_counts_out ${prefix}.ident.csv \\
        ${idents_arg}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    echo '{}' > ${prefix}.xic_config.json
    touch ${prefix}.ident.csv
    """
}
