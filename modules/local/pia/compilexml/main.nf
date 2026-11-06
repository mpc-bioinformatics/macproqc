process PIA_COMPILEXML {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pia:1.5.9--hdfd78af_0':
        'quay.io/biocontainers/pia:1.5.9--hdfd78af_0' }"

    input:
    tuple val(meta), path(identifications)

    output:
    tuple val(meta), path("*.pia.xml"), emit: pia_xml
    tuple val("${task.process}"), val('pia'), eval('pia --version 1>&1 | head -1  | sed "s/.*version //"'), topic: versions, emit: versions_pia

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pia_ram = task.memory.toMega()
    def pia_threads = task.cpus
    """
    pia \\
        ${args} \\
        -Xms2g \\
        -Xmx${pia_ram}m \\
        --threads ${pia_threads} \\
        --compile \\
        -o '${prefix}.pia.xml' \\
        '${identifications}'
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pia.xml
    """
}
