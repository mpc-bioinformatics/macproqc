/*
 * Convert Bruker .d-Folder to mzML using tdf2mzml
 **/
process TDF2MZML {
    tag "${meta.id}"
    label 'process_low'

    container "docker.io/mfreitas/tdf2mzml:0.6.1_noentry"

    input:
    tuple val(meta), path(d_folder)

    output:
    tuple val(meta), path("*.mzML"), emit: spectra
    tuple val("${task.process}"), val('tdf2mzml'), eval('tdf2mzml --version | sed "s/.*tdf2mzml //"'), topic: versions, emit: versions_tdf2mzml

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"


    """
    export MKL_NUM_THREADS=${task.cpus}
    export NUMEXPR_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}

    tdf2mzml ${args} \
            -i ${d_folder} \
            --compression "zlib" \
            -o ${prefix}.mzML
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo ${args}

    touch ${prefix}.mzML
    """
}
