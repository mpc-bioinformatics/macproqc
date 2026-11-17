#!/usr/bin/env nextflow

// 
// Workflows and processes for the export of metrics in mzQC
//

workflow hdf5_to_mzqc {
    take:
        hdf5_files
        output_folder
    
    main:
        mzqc_files = hdf5_to_mzqc_export(hdf5_files.flatten(), output_folder)

    emit:
        mzqc_files.mzqc
}


process hdf5_to_mzqc_export {
    label 'mcquac_image'

    cpus 2
    memory "8.GB"

    publishDir path: { "${output_folder}/mzQC" }, mode: 'copy', pattern: "*.mzqc"

    input:
    path hdf5_files
    val output_folder

    output:
    path "*.mzqc", emit: mzqc

    script:
    """
    hdf5_to_mzqc.py -hdf5 ${hdf5_files} -mzqc_out ${hdf5_files.baseName}.mzqc
    """
}