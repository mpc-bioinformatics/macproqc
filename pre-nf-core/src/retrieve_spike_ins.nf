#!/usr/bin/env nextflow

// 
// Workflows and processees for the extraction of spike in metrics
// 

/*
 * Extracts for each defined spike-in the XIC in the respective m/z and RT region, looks
 * for identifications of the respective sequence and returns the highest peak of the XIC,
 * its RT and the number of identifications.
 *
 * @return spike_in_metrics the metrics of the spike-ins
 */
workflow retrieve_spike_ins_information {
    take:
        raw_files
        psm_mztab_files
        spike_ins_table
        max_parallel_xic_extractors_factor
    
    main:
        // Finally, we generate the input json, retrieve it via trfp and parse back this results into a csv-format
        json_and_identifications = generate_json_and_identifications(psm_mztab_files, spike_ins_table)

        // Map the run basename to the raw files and to their corresponding mzTab files
        runbase_to_raw_and_json = 
            raw_files.map {
                file -> tuple(file.baseName, file)
            }.join(
                json_and_identifications.json.map{
                    file -> tuple(file.name.take(file.name.lastIndexOf('-trfp_input.json')), file)
                },
                by: 0
            )
        
        // branch to Thermo's .raw and Bruker's .d
        runbase_to_raw_and_json.branch {
            thermo: it[1].getExtension() == 'raw'
            bruker: it[1].getExtension() == 'd'
        }.set{ branched_runbase_to_raw_and_json }

        thermo_xics = retrieve_xics_from_thermo_raw_spectra(branched_runbase_to_raw_and_json.thermo, max_parallel_xic_extractors_factor)
        bruker_xics = retrieve_xics_from_bruker_raw_spectra(branched_runbase_to_raw_and_json.bruker, max_parallel_xic_extractors_factor)
        
        xics = thermo_xics.concat(bruker_xics)

        runbase_to_xics_and_identifications = 
            xics.map {
                file -> tuple(file.baseName, file)
            }.join(
                json_and_identifications.identifications.map{
                    file -> tuple(file.name.take(file.name.lastIndexOf('-identifications.csv')), file)
                },
                by: 0
            )

        spike_in_metrics = get_spike_in_metrics(runbase_to_xics_and_identifications, spike_ins_table)

    emit:
        spike_in_metrics
}

/*
 * Creates a json file for the XIC extraction of the spike-ins. For this, the
 * spike-ins file and the identification results are used to refine the retention
 * times for the extraction.
 *
 * @return *-trfp_input.json the JSON for the XIC extraction
 * @return *-identifications.csv a file mapping from the sequences to found identifications
 */
process generate_json_and_identifications {
    label 'mcquac_image'

    cpus 1
    memory '1.GB'

    input:
    path psm_mztab_files
    path spike_ins

    output:
    path "${psm_mztab_files.name.take(psm_mztab_files.name.lastIndexOf('-piaExport-PSM.mzTab'))}-trfp_input.json", emit: json
    path "${psm_mztab_files.name.take(psm_mztab_files.name.lastIndexOf('-piaExport-PSM.mzTab'))}-identifications.csv", emit: identifications

    script:
    """
    create_spike_in_xic_json.py -icsv ${spike_ins} -iidents ${psm_mztab_files} \
            -ojson ${psm_mztab_files.name.take(psm_mztab_files.name.lastIndexOf('-piaExport-PSM.mzTab'))}-trfp_input.json \
            -oidentifications ${psm_mztab_files.name.take(psm_mztab_files.name.lastIndexOf('-piaExport-PSM.mzTab'))}-identifications.csv
    """
}

/**
 * Performs the XIC extraction on Thermo raw files
 *
 * @return baseName.json the extracted XICs in JSON format
 */ 
process retrieve_xics_from_thermo_raw_spectra {
    label 'thermorawfileparser_image'

    cpus {max_parallel_xic_extractors_factor}
    memory '8.GB'

    input:
    tuple val(run_basename), path(raw), path(xic_json)
    val max_parallel_xic_extractors_factor

    output:
    path "${raw.baseName}.json"

    script:
    """
    thermorawfileparser xic -i ${raw} -j ${xic_json}
    """
}

/**
 * Performs the XIC extraction on Bruker .d files
 *
 * @return baseName.json the extracted XICs in JSON format
 */ 
process retrieve_xics_from_bruker_raw_spectra {
    label 'mcquac_image'
    
    cpus {max_parallel_xic_extractors_factor}
    memory '8.GB'

    input:
    tuple val(run_basename), path(raw), path(xic_json)
    val max_parallel_xic_extractors_factor

    output:
    path "${raw.baseName}.json"

    script:
    """
    extract_xic_bruker.py -d_folder ${raw} -in_json ${xic_json} -out_json ${raw.baseName}.json
    """
}

/**
 * Creates the metrics of the spike in extraction, regarding XICs, identification and RT deltas
 *
 * @return baseName-spikeins.csv the metrics of the spike-ins extraction, in CSV
 */ 
process get_spike_in_metrics {
    label 'mcquac_image'

    cpus 1
    memory '4.GB'

    input:
    tuple val(run_basename), path(xic_json), path(identifications)
    path spike_ins_table

    output:
    path "${xic_json.baseName}-spikeins.hdf5"

    script:
    """
    extract_spike_metrics.py -itrfp_json ${xic_json} -iidentifications ${identifications} -ispikeins ${spike_ins_table} -ohdf5 ${xic_json.baseName}-spikeins.hdf5
    """
}
