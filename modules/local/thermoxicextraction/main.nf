/*
 * Extract XICs for spike-in/reference peptides from Thermo raw files, based on a
 * ThermoRawFileParser XIC extraction configuration (see XICEXTRACTIONCONFIG)
 **/
process THERMOXICEXTRACTION {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/thermorawfileparser:1.4.5--h05cac1d_1' :
        'quay.io/biocontainers/thermorawfileparser:1.4.5--h05cac1d_1' }"

    input:
    tuple val(meta), path(raw), path(xic_config)

    output:
    tuple val(meta), path("*.json"), emit: xic
    tuple val("${task.process}"), val('thermorawfileparser'), eval("ThermoRawFileParser.sh --version"), topic: versions, emit: versions_thermoxicextraction

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ThermoRawFileParser.sh xic \\
        ${args} \\
        -i ${raw} \\
        -j ${xic_config}
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo ${args}

    echo '{}' > ${raw.baseName}.json
    """
}
