#!/usr/bin/env nextflow

// 
// Workflows for collecting various metadata of the MS runs from RAW files and mzML files.
// e.g. pump preassure, ion injection time, ...
// 

/**
 * Get metadata headers from Thermo and Bruker raw files, like 
 * - Thermo: Lock Mass, Ion Injection Time, 
 * - Bruker: Vacuum_CurrentFore, Vacuum_CurrentHigh, etc.
 *
 * @return Channel of headers HDF5 files
 */
workflow get_headers {
    take:
        thermo_raw_files
        thermo_memory_limit
        thermo_headers
        bruker_raw_files
        bruker_memory_limit
        bruker_headers_to_parse
        bruker_calibrants

    main:
        thermo_headers = extract_headers_from_thermo_raw_files(thermo_raw_files, thermo_memory_limit, thermo_headers)
        bruker_headers = extract_headers_from_bruker_raw_files(bruker_raw_files, bruker_memory_limit, bruker_headers_to_parse, bruker_calibrants)
        headers = thermo_headers.concat(bruker_headers)

    emit:
        headers
}

/**
 * Get additional data from mzML, e.g. MS1_Density-quartiles, MS2_Density-quartiles, RT-TIC-quartiles, ...
 *
 * @param mzmlfiles Channel of mzML files
 * @return Channel of headers HDF5 files
 */
workflow get_mzml_infos {
    take:
        mzmlfiles
        memory_limit
        base_peak_tic_up_to
        filter_threshold
        report_up_to_charge

    main:
        informations = extract_data_from_mzml(mzmlfiles, memory_limit, base_peak_tic_up_to, filter_threshold, report_up_to_charge)

    emit:
        informations
}

/**
 * Get metadata headers from Thermo and Bruker raw files, like 
 *
 * @param thermo_raw_files Channel of Thermo raw files
 * @return HDF5 with extracted headers
 */
process extract_headers_from_thermo_raw_files {
    label 'mcquac_image'
    
    cpus 1
    memory "${memory_limit}"

    errorStrategy 'ignore'
    stageInMode 'copy'

    input:
    path raw
    val memory_limit
    val thermo_headers

    output:
    path "${raw.baseName}-custom_headers.hdf5"

    script:
    """
    # Pythonnet sometimes fails to exit and throws a mono error
    extract_thermo_headers.py -raw ${raw} -out_hdf5 ${raw.baseName}-custom_headers.hdf5 ${thermo_headers} || true

    # Fail Check if no content was written
    if ! [ -s "${raw.baseName}-custom_headers.hdf5" ];then
        rm ${raw.baseName}-custom_headers.hdf5
    fi
    """
}

/**
 * Get metadata headers from Thermo and Bruker raw files, like 
 *
 * @param bruker_raw_files Channel of Bruker .d-folders
 * @return HDF5 with extracted headers
 */
process extract_headers_from_bruker_raw_files {
    label 'mcquac_image'

    cpus 1
    memory "${memory_limit}"

    errorStrategy 'ignore'

    input:
    path raw
    val memory_limit
    val bruker_headers_to_parse
    val bruker_calibrants

    output:
    path "${raw.baseName}-custom_headers.hdf5"

    script:
    """
    extract_bruker_headers.py -d_folder ${raw} -out_hdf5 ${raw.baseName}-custom_headers.hdf5 ${bruker_headers_to_parse} ${bruker_calibrants}
    """
}

/**
 * Get additional data from mzML, e.g. MS1_Density-quartiles, MS2_Density-quartiles, RT-TIC-quartiles, ...
 *
 * @param mzml Channel of mzML files
 * @return HDF5 file with the extracted data
 */
process extract_data_from_mzml {
    label 'mcquac_image'

    cpus 1
    memory "${memory_limit}"

    input:
    path mzml
    val memory_limit
    val base_peak_tic_up_to
    val filter_threshold
    val report_up_to_charge

    output:
    path "${mzml.baseName}-mzml_info.hdf5"

    script:
    """
    extract_data_from_mzml.py -mzml ${mzml} -out_hdf5 ${mzml.baseName}-mzml_info.hdf5 -base_peak_tic_up_to ${base_peak_tic_up_to} -filter_threshold ${filter_threshold} -report_up_to_charge ${report_up_to_charge}
    """
}
