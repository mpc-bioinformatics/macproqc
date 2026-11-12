/*
 * Compute DIA quality metrics (per run and per isolation window) using OpenMS DIAuditor
 *
 * DIAuditor is available from OpenMS 3.6.0 on; there is no bioconda build yet, hence no conda directive.
 **/
process DIAUDITOR {
    tag "${meta.id}"
    label 'process_low'

    container "ghcr.io/openms/openms-tools-thirdparty:v3.6.0"

    input:
    tuple val(meta), path(spectra)   // one or more DIA runs (mzML, Thermo .raw or Bruker .d); one table row per run

    output:
    tuple val(meta), path("*.tsv"), emit: report
    tuple val(meta), path("*_windows.tsv"), emit: windows
    tuple val(meta), path("*.mzQC"), emit: mzqc
    tuple val("${task.process}"), val('openms'), eval("FileInfo --help 2>&1 | sed -nE 's/^Version: ([0-9.]+).*/\\1/p'"), topic: versions, emit: versions_openms

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    DIAuditor \\
        -in ${spectra} \\
        -out ${prefix}.tsv \\
        -out_windows ${prefix}_windows.tsv \\
        -out_mzqc ${prefix}.mzQC \\
        -threads ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.tsv
    touch ${prefix}_windows.tsv
    touch ${prefix}.mzQC
    """
}
