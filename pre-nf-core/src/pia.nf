#!/usr/bin/env nextflow

// 
// Workflows and processes for protein inference using PIA
// 

// Note for CPUs: PIA needs/takes very shortly all available processors, but rather idles on them later.
// On the memory-side, it uses a lot the more PSMs are found...
// Use this information to adjust number of parallel threads.

/*
 * Performs PIA's FDR and protein iference on the given files and levels.
 *
 * @return return_files tuples containing the PSM, peptide and protein level results each (may be empty, if level was not returned)
 */
workflow pia_analysis {
    take:
        identifications
        do_psm_export
        do_peptide_export
        do_protein_export
        fdr_filter
        pia_threads
        pia_gb_ram
        pia_fdr_threshold
        pia_prefilter_threshold

    main:
        basename_to_ids = identifications.map { it -> 
            tuple(
                it.baseName,
                it
            )
        }
        basename_to_pia_xmls = compile_pia_xmls(basename_to_ids, pia_threads, pia_gb_ram)

        if (pia_prefilter_threshold > 0) {
            // perform the pre-filtering
            basename_to_pia_xmls = perform_pia_prefiltering(basename_to_pia_xmls, pia_prefilter_threshold, pia_threads, pia_gb_ram)
        }

        analysis_json = prepare_analysis_json(do_psm_export, do_peptide_export, do_protein_export, fdr_filter,
                true,      // remove_decoys
                pia_fdr_threshold
        )
        pia_all_report_files = pia_run_analysis(basename_to_pia_xmls, analysis_json, pia_threads, pia_gb_ram, "mzTab")

        // TODO: remove this ugly hack! implement basenames and results further downstream correctly!
        pia_all_report_files = pia_rename_files(pia_all_report_files.psms, pia_all_report_files.peptides, pia_all_report_files.proteins)

        return_files = pia_all_report_files.psms.collect()
            .concat(pia_all_report_files.peptides.collect())
            .concat(pia_all_report_files.proteins.collect())
            .toList()
            .transpose()
    
    emit:
        return_files
}


/*
 * Performs PIA's FDR and protein iference on the given files and all levels.
 * 
 * @return return_files tuples containing the PSM, peptide and protein level results each (may be empty, if level was not returned)
 */
workflow pia_analysis_full {
    take:
        identifications
        pia_threads
        pia_gb_ram
        pia_fdr_threshold
        pia_prefilter_threshold

    main:
        pia_report_files = pia_analysis(
            identifications,
            true,
            true,
            true,
            true,
            pia_threads,
            pia_gb_ram,
            pia_fdr_threshold,
            pia_prefilter_threshold
        )
    
    emit:
        pia_report_files
}

/*
 * Performs PIA's FDR and protein iference on the given files and PSM level only.
 *
 * @return psm_file the results on PSM level
 */
workflow pia_analysis_psm_only {
    take:
        identifications
        fdr_filter
        pia_threads
        pia_gb_ram
        pia_fdr_threshold
    
    main:
        pia_report_files = pia_analysis(
            identifications,
            true,
            false,
            false,
            fdr_filter,
            pia_threads,
            pia_gb_ram,
            pia_fdr_threshold,
            0                           // 0 means no pre-filtering here
        )

    emit:
        pia_report_files
            .toList()
            .transpose()
            .first()
            .flatten()
}

/**
 * Extracts the QC metrics from the PIA results and writes them into a HFD5.
 * If peptides_table is true, the peptide tables will be saved to the results.
 *
 * @return extracted hdf5 metrics written into the HFD5 file
 */
workflow pia_extract_metrics {
    take:
        pia_results
        peptides_table
        output_folder

    main:
        extract_csv = pia_extract_csv(pia_results)

        if (peptides_table) {
            save_peptide_files(pia_results, output_folder)
        }   
    
    emit:
        extract_csv
}


workflow perform_pia_prefiltering {
    take:
        basename_to_pia_xmls
        pia_prefilter_threshold
        pia_threads
        pia_gb_ram

    main:
        // perform the pre-filtering
        prefilter_analysis_json = prepare_analysis_json(
            true,  // psm_export
            false, // peptide_export
            false, // protein_export
            true,  // fdr_filter
            false, // remove_decoys
            pia_prefilter_threshold
        )
        pia_prefiltered = pia_run_analysis(basename_to_pia_xmls, prefilter_analysis_json, pia_threads, pia_gb_ram, "mzid")
        prefiltered_pia_xmls = compile_pia_xmls(pia_prefiltered.psms, pia_threads, pia_gb_ram)

    emit:
        prefiltered_pia_xmls
}


/**
 * Compiles the input files into PIA intermediate files
 *
 * @return pia.xml the compilation as pia.xml
 */
process compile_pia_xmls {
    label 'pia_image'
    
    cpus {pia_threads}
    memory "${pia_gb_ram}.GB"

    input:
    tuple val(id_basename), path(identifications)      // mapping from the basename to the identifications
    val pia_threads
    val pia_gb_ram

    output:
    tuple val(id_basename), path("${id_basename}.pia.xml")

    script:
    """
    pia -Xms2g -Xmx${pia_gb_ram}g --threads ${pia_threads} --compile -o '${id_basename}.pia.xml' '${identifications}'
    """
}

/**
 * Creates a PIA json analysis file with given export parameters
 *
 * @return an pia_analysis.json with defined parameters for the QC
 */
process prepare_analysis_json {
    label 'pia_image'

    cpus 1
    memory "4.GB"

    input:
    val psm_export
    val peptide_export
    val protein_export
    val fdr_filter
    val remove_decoys       // remove decoys works only together with fdr_filter == true
    val fdr_threshold

    output:
    path "pia_analysis.json"

    script:
    """
    pia --example > pia_analysis.json
    
    # delete the first row of the file, which contains an explanatory comment and does not start with the json
    sed -i 1d pia_analysis.json

    sed -i 's;"createPSMsets": .*,;"createPSMsets": false,;g' pia_analysis.json
    sed -i 's;"psmLevelFileID": .*,;"psmLevelFileID": 1,;g' pia_analysis.json
    sed -i 's;"topIdentifications": .*,;"topIdentifications": 1,;g' pia_analysis.json
    sed -i 's;"calculateCombinedFDRScore": .*,;"calculateCombinedFDRScore": false,;g' pia_analysis.json
    if [ ${psm_export} = true ];
    then
      sed -i 's;"psmExportFile":.*,;"psmExportFile": "piaExport-PSMs.mzTab",;g' pia_analysis.json
    else
      sed -i 's;"psmExportFile":.*;;g' pia_analysis.json
    fi

    if [ ${peptide_export} = true ];
    then
      sed -i 's;"inferePeptides":.*,;"inferePeptides": true,;g' pia_analysis.json
      sed -i 's;"peptideExportFile":.*,;"peptideExportFile": "piaExport-peptides.csv",;g' pia_analysis.json
    else
      sed -i 's;"inferePeptides":.*,;"inferePeptides": false,;g' pia_analysis.json
      sed -i 's;"peptideExportFile":.*,;;g' pia_analysis.json
    fi
    sed -i 's;"peptideLevelFileID":.*,;"peptideLevelFileID": 1,;g' pia_analysis.json

    if [ ${protein_export} = true ];
    then
      sed -i 's;"infereProteins":.*,;"infereProteins": true,;g' pia_analysis.json
      sed -i 's;"inferenceMethod":.*,;"inferenceMethod": "inference_spectrum_extractor",;g' pia_analysis.json
      sed -i '/inferenceFilters/{n;s/.*/    "psm_score_filter_psm_fdr_score <= 0.01"/g;}' pia_analysis.json
      sed -i 's;"scoringBaseScore":.*,;"scoringBaseScore": "psm_fdr_score",;g' pia_analysis.json
      sed -i 's;"scoringPSMs":.*,;"scoringPSMs": "best",;g' pia_analysis.json
      sed -i 's;"proteinExportFile":.*,;"proteinExportFile": "piaExport-proteins.mzTab",;g' pia_analysis.json
      sed -i 's;"proteinExportWithPSMs":.*,;"proteinExportWithPSMs": true,;g' pia_analysis.json
    else
      sed -i 's;"infereProteins":.*,;"infereProteins": false,;g' pia_analysis.json
      sed -i 's;"proteinExportFile":.*,;;g' pia_analysis.json
    fi

    if [ ${fdr_filter} = true ];
    then
      if [ ${remove_decoys} = true ];
      then
        sed -i '/psmFilters/{n;s/.*/    "psm_score_filter_psm_fdr_score <= ${fdr_threshold}",\\n    "psm_accessions_filter !regex_only DECOY_.*"/g;}' pia_analysis.json
        sed -i '/peptideFilters/{n;s/.*/    "psm_score_filter_psm_fdr_score <= ${fdr_threshold}",\\n    "peptide_accessions_filter !regex_only DECOY_.*"/g;}' pia_analysis.json
        sed -i '/proteinFilters/{n;s/.*/    "protein_q_value_filter <= ${fdr_threshold}",\\n    "protein_accessions_filter !regex_only DECOY_.*"/g;}' pia_analysis.json
      else
        sed -i '/psmFilters/{n;s/.*/    "psm_score_filter_psm_fdr_score <= ${fdr_threshold}"/g;}' pia_analysis.json
        sed -i '/peptideFilters/{n;s/.*/    "psm_score_filter_psm_fdr_score <= ${fdr_threshold}"/g;}' pia_analysis.json
        sed -i '/proteinFilters/{n;s/.*/    "protein_q_value_filter <= ${fdr_threshold}"/g;}' pia_analysis.json
      fi
    else
      sed -i '/psmFilters/{n;s/.*//g;}' pia_analysis.json
      sed -i '/peptideFilters/{n;s/.*//g;}' pia_analysis.json
      sed -i '/proteinFilters/{n;s/.*//g;}' pia_analysis.json
    fi
    """
}


/*
 * Runs a PIA analysis with a parameter file
 * The command line allows you to execute an analysis via prior defined analysis in JSON format. 
 * Additionally to the json file, the prior compiled intermediate file must be given.
 *
 * @return exports for the PSM, peptide and protein level (might be emmpty files, if teh analysis is given like it)
 *
 */
process pia_run_analysis {
    label 'pia_image'

    cpus {pia_threads}
    memory "${pia_gb_ram}.GB"

    input:
    tuple val(search_basename), path(pia_xml)
    path pia_analysis_file
    val pia_threads
    val pia_gb_ram
    val psm_export_format               // needs to be set correctly for pre-processing (mzid) and final (mzTab) 

    output:
    tuple val(search_basename), path("psms.${psm_export_format}"), emit: psms
    tuple val(search_basename), path("peptides.csv"), emit: peptides
    tuple val(search_basename), path("proteins.mzTab"), emit: proteins

    script:
    """
    jsonfile="${search_basename}-analysis.json"

    cp ${pia_analysis_file} \${jsonfile}

    touch psms.${psm_export_format}
    touch peptides.csv
    touch proteins.mzTab

    sed -i 's;"psmExportFile": .*;"psmExportFile": "psms.${psm_export_format}",;g' \${jsonfile}
    sed -i 's;"peptideExportFile": .*;"peptideExportFile": "peptides.csv",;g' \${jsonfile}
    sed -i 's;"proteinExportFile": .*;"proteinExportFile": "proteins.mzTab",;g' \${jsonfile}

    pia -Xms2g -Xmx${pia_gb_ram}g --threads ${pia_threads} \${jsonfile} '${pia_xml}'
    """
}


process pia_extract_csv {
    label 'mcquac_image'

    cpus 1
    memory "8.GB"

    input:
    tuple path(psm_results), path(peptide_results), path(protein_results)

    output:
    path "${psm_results.name.take(psm_results.name.lastIndexOf('-piaExport-PSM.mzTab'))}-pia_extraction.hdf5"

    script:
    """
    outfile="${psm_results.name.take(psm_results.name.lastIndexOf('-piaExport-PSM.mzTab'))}-pia_extraction.hdf5"
    extract_from_pia_output.py --pia_PSMs ${psm_results} --pia_peptides ${peptide_results} --pia_proteins ${protein_results} --out_hdf5 \${outfile}
    """
}


process save_peptide_files {
    label 'mcquac_image'

    cpus 1
    memory "1.GB"

    publishDir path: { "${output_folder}" }, mode: 'move', overwrite: true

    input:
    tuple path(psm_results), path(peptide_results), path(protein_results)
    val output_folder

    output:
    path "peptides/${peptide_results.name}"

    script:
    """
    # create a symlink to save
    mkdir peptides
    cp ${peptide_results} peptides/${peptide_results.name}
    """
}


//
// TODO: GET RID OF THIS PROCESS and implement the file and basenames correctly further downstream
//
process pia_rename_files {
    label 'mcquac_image'

    cpus 1
    memory "4.GB"

    input:
    tuple val(psms_basename), path(psms_result)
    tuple val(peptides_basename), path(peptides_result)
    tuple val(proteins_basename), path(proteins_result)

    output:
    path "${psms_basename}-piaExport-PSM.mzTab", emit: psms
    path "${peptides_basename}-piaExport-peptides.csv", emit: peptides
    path "${proteins_basename}-piaExport-proteins.mzTab", emit: proteins

    script:
    """
    cp ${psms_result} ${psms_basename}-piaExport-PSM.mzTab
    cp ${peptides_result} ${peptides_basename}-piaExport-peptides.csv
    cp ${proteins_result} ${proteins_basename}-piaExport-proteins.mzTab
    """
}
