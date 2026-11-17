#!/usr/bin/env nextflow

// 
// Workflows and processes for handling process exit informations.
// Makes sure the user knows if a file was processed successfully or not.
// 


/*
 * Writes out processing information
 * 
 * @param raw_files the raw files of the complete project
 * @param processed_basenames the basenames of the successfully processed files
 *
 * @return collected_information a file containing human readable information about whether the files were processed successfully
 */
workflow output_processing_success {
    take:
    raw_files
    processed_basenames
    output_dir

    main:
    raw_map = raw_files.map{file -> tuple(file.baseName, file)}
	processed_ok = processed_basenames.map{name -> tuple(name, "successful")}

    info_files = write_processing_information(raw_map.join(processed_ok, remainder: true))
    
    collected_information = info_files.collectFile(
            storeDir: "${output_dir}/",
            name:  'processing_info.txt',
            newLine: true
        )
    
    emit:
    collected_information
}


/**
 * Writes for each raw file whether it was processed successfully or not.  
 *
 * @param tuple of val(file_base), path(raw_file), val(information), while information is "successful" or null
 *
 * @return a file containing human readable info about the processing state
 */
process write_processing_information {
    label 'mcquac_image'

    cpus 1
    memory '1.GB'

	input:
	tuple val(file_base), path(raw_file), val(information)

    output:
    path "${file_base}_processing_info.txt"

    script:
    """
    if [ "${information}" = "null" ]; then
        msg="ERROR processing"
    else
        msg="successfully processed"
    fi
    
    echo "\${msg} ${file_base} (${raw_file})" > ${file_base}_processing_info.txt
    """
}