
process QCVISUALIZATION {
    tag "$meta.id"
    label 'process_medium'

    container 'ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-74e33ed'

    input:
    tuple val(meta), path(hdf5_files)   // tuple: meta information and list of hdf5 files containing the metrics to visualize
    val output_table_type               // val: output table type for the QC visualization (tsv, csv, xlsx; default: tsv)
    val figure_format                   // val: figure format for the QC visualization (plotly, html, or both as 'html,plotly'; default: plotly)
    path spike_ins_table                // path: spike-ins reference table as a csv file, or [] if no spike-ins are used

    output:
    tuple val(meta), path("*.json"), emit: jsons, optional: true
    tuple val(meta), path("*.html"), emit: htmls, optional: true
    tuple val(meta), path("*.{csv,tsv,xlsx}"), emit: output_tables, optional: false
    tuple val(meta), path("fig13_MS1_map"), emit: ms1_maps, optional: false
    tuple val(meta), path("fig16_additional_headers"), emit: additional_plots, optional: true
    tuple val(meta), path("fig17_BRUKER_calibrants"), emit: bruker_calibrants, optional: true
    tuple val("${task.process}"), val('macproqc_helpers'), eval('macproqc-helpers --version | sed "s;^macproqc-helpers ;;"'), topic: versions, emit: versions_macproqc_helpers

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''   // free-form, optional extra CLI arguments
    def spike_ins_tab = spike_ins_table ? "-spike_ins_table ${spike_ins_table}" : ''
    hdf5_files_list = hdf5_files.join(' ')
    // additional arguments can be passed with the args variable:
    // -group GROUP                                 List of the experimental group (comma-separated).
    // -RT_unit RT_UNIT                             Unit of the retention time, either sec for seconds or min for minutes.
    // -output_column_order OUTPUT_COLUMN_ORDER     Order of columns in the output table
    // -spikein_columns SPIKEIN_COLUMNS             Columns of the spike-in dataframes that should end up in the result table
    // -height_barplots HEIGHT_BARPLOTS             Height of the barplots in pixels
    // -width_barplots WIDTH_BARPLOTS               Width of the barplots in pixels
    // -height_pca HEIGHT_PCA                       Height of the PCA plots in pixels
    // -width_pca WIDTH_PCA                         Width of the PCA plots in pixels
    // -height_ionmaps HEIGHT_IONMAPS               Height of the ionmaps in inches
    // -width_ionmaps WIDTH_IONMAPS                 Width of the ionmaps in inches
    // -metric_registry_file METRIC_REGISTRY_FILE   Path to a JSON file storing the metric-type registry
    """
    python -m macproqc_helpers visualize \\
        -hdf5_files ${hdf5_files_list} \\
        -output "." \\
        -output_table_type ${output_table_type} \\
        -figure_format ${figure_format} \\
        ${spike_ins_tab} \\
        -log_file \\
        ${args}
    """


    stub:
    """
    touch Figure1.json
    touch Outputtable.csv
    mkdir fig13_MS1_map
    mkdir fig16_additional_headers
    mkdir fig17_BRUKER_calibrants
    """
}
