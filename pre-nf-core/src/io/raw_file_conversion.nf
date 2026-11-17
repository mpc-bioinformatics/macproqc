#!/usr/bin/env nextflow

//
// Workflows and processes file conversions, e.g. from vendor file formats to open standards like mzML 
//

/**
 * Convert raw files (Thermo Fisher .raw-files and Bruker tdf-files) to mzML files
 * @params thermo_raw_files Thermo Fisher .raw-files
 * @params bruker_raw_folders Bruker .d-folders
 */
workflow convert_raws_to_mzml {
    take:
        thermo_raw_files
        bruker_raw_folders
        thermo_conversion_mem
        bruker_conversion_cpu
        bruker_conversion_mem
        
    main:
        thermo_mzmls = convert_thermo_raw_files(thermo_raw_files, thermo_conversion_mem)
        bruker_mzmls = convert_bruker_raw_folders(bruker_raw_folders, bruker_conversion_cpu, bruker_conversion_mem)
        mzmls = thermo_mzmls.concat(bruker_mzmls)
    emit:
        mzmls
}

/**
 * Convert raw file (Thermo Fisher .raw-files) to mzML files
 * @params thermo_raw_files A list of Thermo Fisher .raw-files
 *
 * @return mzML files
 */
process convert_thermo_raw_files {
    label 'thermorawfileparser_image'
    errorStrategy 'ignore'

    // Thermo Raw File parser is currently limited to 2 CPUs, see:
    // * https://github.com/compomics/ThermoRawFileParser/issues/23
    // * https://github.com/compomics/ThermoRawFileParser/issues/95
    cpus 2
    memory "${thermo_conversion_mem}"

    input:
    path raw_file
    val thermo_conversion_mem

    output:
    path "${raw_file.baseName}.mzML"

    script:
    """
    thermorawfileparser --format=2 --output_file=${raw_file.baseName}.mzML --input=${raw_file}
    """
}

/**
 * Convert raw files (Bruker .d-folder) to mzML files
 * @params raw_folder Bruker .d-folder
 *
 * @return mzML file
 */
process convert_bruker_raw_folders {
    label 'tdf2mzml_image'
    errorStrategy 'ignore'
    
    cpus "${bruker_conversion_cpu}"
    memory "${bruker_conversion_mem}"

    input:
    path raw_folder
    val bruker_conversion_cpu
    val bruker_conversion_mem
    
    output:
    path "${raw_folder.baseName}.mzML"

    script:
    """
    export MKL_NUM_THREADS=${bruker_conversion_cpu}
    export NUMEXPR_NUM_THREADS=${bruker_conversion_cpu}
    export OMP_NUM_THREADS=${bruker_conversion_cpu}

    tdf2mzml -i ${raw_folder} --compression "none" --ms1_type centroid -o ${raw_folder.baseName}.mzML
    """
}
