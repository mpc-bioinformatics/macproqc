/*
 * Create a ThermoRawFileParser (also alphatims-compatible) XIC extraction config
 * from a spike-ins CSV and (optionally) a PSM mzTab file to refine the RT windows
 **/
process XICEXTRACTIONCONFIG {
    tag "${meta.id}"
    label 'process_single'

    container "ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-74e33ed"

    input:
    tuple val(meta), path(psm_mztab_file)   // psm_mztab_file: PSM mzTab, or [] to use only the RT windows of the spike-ins table
    path spike_ins

    output:
    tuple val(meta), path("*.json"), path("*.ident.csv"), emit: config
    tuple val("${task.process}"), val('macproqc_helpers'), eval('macproqc-helpers --version | sed "s;^macproqc-helpers ;;"'), topic: versions, emit: versions_macproqc_helpers

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def idents_arg = psm_mztab_file ? "-iidents ${psm_mztab_file}" : ''

    """
    python -m macproqc_helpers create_thermorawfileparser_xic_config \\
        ${args} \\
        -icsv ${spike_ins} \\
        ${idents_arg} \\
        -ojson ${prefix}.json \\
        -oidentifications ${prefix}.ident.csv
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    echo '{}' > ${prefix}.json
    touch ${prefix}.ident.csv
    """
}
