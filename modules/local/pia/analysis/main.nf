process PIA_ANALYSIS {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pia:1.5.9--hdfd78af_0':
        'quay.io/biocontainers/pia:1.5.9--hdfd78af_0' }"

    input:
    tuple val(meta), path(pia_xml), path(pia_json)

    output:
    tuple val(meta), path("*.piaExport-PSMs.{mzTab,mzid}"), optional: true, emit: psms
    tuple val(meta), path("*.piaExport-peptides.csv"), optional: true, emit: peptides
    tuple val(meta), path("*.piaExport-proteins.mzTab"), optional: true, emit: proteins
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
        '${pia_json}' \\
        '${pia_xml}'

    psms_file=(piaExport-PSMs.*)
    if [ -e "\${psms_file[0]}" ]; then
        mv "\${psms_file[0]}" "${prefix}.\${psms_file[0]}"
    fi

    if [ -f piaExport-peptides.csv ]; then
        mv piaExport-peptides.csv ${prefix}.piaExport-peptides.csv
    fi

    if [ -f piaExport-proteins.mzTab ]; then
        mv piaExport-proteins.mzTab ${prefix}.piaExport-proteins.mzTab
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.piaExport-PSMs.mzTab
    touch ${prefix}.piaExport-peptides.csv
    touch ${prefix}.piaExport-proteins.mzTab
    """
}
