//
// Subworkflow with functionality specific to the nf-core/macproqc pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet (only used, if the metrics extraction is performed)
    input_hdf5_dir    //  string: Path to a directory with metric HDF5 files (only used, if the metrics extraction is skipped)
    fasta             //  string: Path to the FASTA database (only used, if the identification is performed)
    skip_metrics_extraction // boolean: Skip the metrics extraction (requires input_hdf5_dir)
    skip_identification     // boolean: Skip the identification during the metrics extraction
    skip_mzqc_export        // boolean: Skip the mzQC export
    skip_visualization      // boolean: Skip the visualization
    spike_ins_table   //  string: Path to the spike-ins table (only used, if the metrics extraction or visualization is performed)
    help             // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/macproqc ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/macproqc/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --fasta database.fasta --outdir <OUTDIR>\n" +
        "  or, to only export/visualize the HDF5 files of an earlier run:\n" +
        "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --skip_metrics_extraction --input_hdf5_dir <OUTDIR>/metrics/hdf5 --outdir <NEW_OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        false
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Check the combination of input parameters and pipeline steps
    //
    validateInputParameters(input, input_hdf5_dir, fasta, skip_metrics_extraction, skip_identification, skip_mzqc_export, skip_visualization, spike_ins_table)

    //
    // Validate the spike-ins table (only if it is used by the metrics extraction or visualization)
    //
    if (spike_ins_table && !(skip_metrics_extraction && skip_visualization)) {
        samplesheetToList(spike_ins_table, "${projectDir}/assets/schema_spike_ins.json")
    }

    //
    // Create channel from input file provided through params.input (for the metrics extraction)
    //
    def ch_samplesheet = channel.empty()
    if (!skip_metrics_extraction) {
        ch_samplesheet = channel
            .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
            .map { row ->
                def meta          = row[0]
                def spectrum_file = row[1]

                return [meta.id, meta, spectrum_file]
            }
            .groupTuple()
            .map { samplesheet -> validateInputSamplesheet(samplesheet) }
    }

    //
    // Create channel from the HDF5 files in params.input_hdf5_dir (if the metrics extraction is skipped)
    // the run id is restored from the file name, as COMBINEHDF5 writes "<id>.hdf5"
    //
    def ch_input_hdf5 = channel.empty()
    if (skip_metrics_extraction) {
        ch_input_hdf5 = channel
            .fromPath("${input_hdf5_dir}/*.hdf5", checkIfExists: true)
            .map { hdf5 -> [ [ id: hdf5.baseName ], hdf5 ] }
    }

    emit:
    samplesheet = ch_samplesheet
    hdf5        = ch_input_hdf5     // channel: [ val(meta), path(hdf5) ]
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                []
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Check the combination of input parameters and enabled pipeline steps
// (--input and --input_hdf5_dir are mutually exclusive, depending on the metrics extraction)
//
def validateInputParameters(input, input_hdf5_dir, fasta, skip_metrics_extraction, skip_identification, skip_mzqc_export, skip_visualization, spike_ins_table) {
    if (skip_metrics_extraction && skip_mzqc_export && skip_visualization) {
        error("All pipeline steps are skipped (--skip_metrics_extraction, --skip_mzqc_export and --skip_visualization), nothing to do.")
    }

    if (!skip_metrics_extraction) {
        if (!input) {
            error("Please provide a samplesheet with --input (or skip the metrics extraction with --skip_metrics_extraction and use --input_hdf5_dir).")
        }
        if (input_hdf5_dir) {
            error("--input_hdf5_dir can only be used together with --skip_metrics_extraction (and without --input).")
        }
        if (!skip_identification && !fasta) {
            error("Please provide a FASTA database with --fasta (or skip the identification with --skip_identification).")
        }
    } else {
        if (!input_hdf5_dir) {
            error("--skip_metrics_extraction requires the directory of earlier extracted HDF5 files given with --input_hdf5_dir.")
        }
        if (input) {
            error("--input cannot be used together with --skip_metrics_extraction, use --input_hdf5_dir for the HDF5 files instead.")
        }
        if (skip_identification) {
            log.warn("--skip_identification has no effect, as the metrics extraction is skipped.")
        }
    }

    if (fasta && (skip_metrics_extraction || skip_identification)) {
        log.warn("--fasta is ignored, as the identification is not performed.")
    }

    if (spike_ins_table && skip_metrics_extraction && skip_visualization) {
        log.warn("--spike_ins_table is ignored, as neither the metrics extraction nor the visualization is performed.")
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, spectra_files) = input[1..2]

    if (spectra_files.size() != 1) {
        error("Please check input samplesheet -> Expected exactly one 'spectrum_file' per sample row after grouping, got ${spectra_files.size()} for: ${metas[0].id}")
    }

    return [metas[0], spectra_files[0]]
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
