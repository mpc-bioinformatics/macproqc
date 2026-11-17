process COMETCONFIG {
    tag "$meta.id"
    label 'process_single'

    publishDir path: { "${params.outdir}/comet" }, mode: params.publish_dir_mode

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.14'
        : 'biocontainers/python:3.14'}"

    input:
    tuple val(meta), path(comet_config_template)

    output:
    tuple val(meta), path("${params_out_file}"), emit: params
    path "versions.yml"                        , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Supported ext.args (all optional; units and ranges match the Comet params file):
    //   -peptide_mass_tolerance_upper <float>  (default: 5.0)
    //   -peptide_mass_tolerance_lower <float>  (default: -5.0)  — usually negative
    //   -peptide_mass_units           <int>    (default: 2)      0=amu, 1=mmu, 2=ppm
    //   -isotope_error                <int>    (default: 2)      0=off, 1=0/1, 2=0/1/2, 3=0/1/2/3, 4=-1/0/1/2/3, 5=-1/0/1
    //   -fragment_bin_tol             <float>  (default: 0.02)
    //   -fragment_bin_offset          <float>  (default: 0.0)    0.0–1.0
    //   -theoretical_fragment_ions    <int>    (default: 0)      0=flanking peaks, 1=M peak only
    //   -psms_per_spectrum            <int>    (default: 5)      maps to num_output_lines
    //   -variable_modifications       <string> (default: "")     semicolon-separated Comet variable_modXX entries
    //   -static_modifications         <string> (default: "")     semicolon-separated Comet add_XX entries
    def prefix = task.ext.prefix ?: "${meta.id}"
    params_out_file = "${prefix}.comet.params"
    template 'adjust_comet_params.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    params_out_file = "${prefix}.comet.params"
    """
    touch ${params_out_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
