process IDMETRICSEXTRACTION {
    tag "$meta.id"
    label 'process_single'

    container 'ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-55b549b'

    input:
    tuple val(meta), path(psms), path(peptides), path(proteins)

    output:
    tuple val(meta), path("${prefix}.idmetrics.hdf5"), emit: hdf5
    tuple val("${task.process}"), val('macproqc_helpers'), eval('macproqc-helpers --version | sed "s;^macproqc-helpers ;;"'), topic: versions, emit: versions_macproqc_helpers

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // missing PIA results ([]) are replaced by empty files, which ultimately returns 0 identifications
    def psms_file     = psms     ?: 'no_results.piaExport-PSMs.mzTab'
    def peptides_file = peptides ?: 'no_results.piaExport-peptides.csv'
    def proteins_file = proteins ?: 'no_results.piaExport-proteins.mzTab'
    def empty_files   = [ psms ? null : psms_file, peptides ? null : peptides_file, proteins ? null : proteins_file ].findAll()
    """
    # create empty files for missing PIA results
    ${empty_files ? "touch ${empty_files.join(' ')}" : ''}

    python -m macproqc_helpers collect-pia-metrics \\
        --pia_PSMs ${psms_file} \\
        --pia_peptides ${peptides_file} \\
        --pia_proteins ${proteins_file} \\
        --out_hdf5 ${prefix}.idmetrics.hdf5 \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.idmetrics.hdf5
    """
}
