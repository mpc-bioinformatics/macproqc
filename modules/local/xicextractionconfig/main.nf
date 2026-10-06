/*
 * Create a ThermoRawFileParser (also alphatims-compatible) XIC extraction config
 * from a PSM mzTab file and a spike-ins CSV
 **/
process XICEXTRACTIONCONFIG {
    tag "${meta.id}"
    label 'process_single'

    container "ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-a9c31aa"

    input:
    tuple val(meta), path(psm_mztab_file)
    path spike_ins

    output:
    tuple val(meta), path("*.json"), path("*.ident.csv"), emit: config
    tuple val("${task.process}"), val('macproqc_helpers'), eval('macproqc-helpers --version | sed "s;^macproqc-helpers ;;"'), topic: versions, emit: versions_macproqc_helpers

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    python -m macproqc_helpers create_thermorawfileparser_xic_config \\
        ${args} \\
        -icsv ${spike_ins} \\
        -iidents ${psm_mztab_file} \\
        -ojson ${psm_mztab_file}.json \\
        -oidentifications ${psm_mztab_file}.ident.csv
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
