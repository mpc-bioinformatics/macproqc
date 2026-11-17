#!/usr/bin/env nextflow

// 
// Workflows and processes for peptide identification using Comet
// 


/*
 * Identifies peptides in MS/MS spectra using Comet
 *
 * @return mzID containing the identification results
 */
workflow identification_with_comet {
    take:
        mzmls
        fasta_file
        generate_decoys
        decoy_method
        fasta_output_folder
        store_decoy_fasta
        comet_cpu
        comet_mem
        peptide_mass_tolerance_upper
        peptide_mass_tolerance_lower
        peptide_mass_units
        isotope_error
        fragment_bin_tol
        fragment_bin_offset
        theoretical_fragment_ions
        static_modifications

    main:
        comet_params_file = create_fresh_comet_params()
        adj_comet_params = adjust_comet_params(
            comet_params_file,
            comet_cpu,
            peptide_mass_tolerance_upper,
            peptide_mass_tolerance_lower,
            peptide_mass_units,
            isotope_error,
            fragment_bin_tol,
            fragment_bin_offset,
            theoretical_fragment_ions,
            static_modifications
        )

        if (generate_decoys) {
            fasta_file = generate_decoy_database(
                fasta_file,
                decoy_method,
                store_decoy_fasta,
                fasta_output_folder
            )
        }

        id_results = comet_search(
            mzmls,
            fasta_file,
            adj_comet_params,
            comet_cpu,
            comet_mem
        )
    
    emit:
        mzids = id_results.mzids
}

/**
 * Creates a fresh comet params file
 */
process create_fresh_comet_params {
    label 'comet_image'

    cpus 1
    memory "1.GB"

    output:
    path "comet.params"

    script:
    """
    # create a new comet params file
    comet -p
    mv comet.params.new comet.params
    """
}

/**
 * Adjust comet.params to have the correct output files and limited allowed threads.
 * Also set all parameters.
 *
 * @return adjusted comet parameter file
 */
process adjust_comet_params {
    label 'mcquac_image'

    cpus 1
    memory "1.GB"

    input:
    path comet_params_file
    val comet_cpu
    val peptide_mass_tolerance_upper
    val peptide_mass_tolerance_lower
    val peptide_mass_units
    val isotope_error
    val fragment_bin_tol
    val fragment_bin_offset
    val theoretical_fragment_ions
    val static_modifications

    output:
    path "adjusted.comet.params"

    script:
    """
    # set the number of threads
    sed -i 's/^num_threads.*/num_threads = ${comet_cpu}/' ${comet_params_file}

    adjust_comet_params.py -comet_params ${comet_params_file} -params_out adjusted.comet.params \
        -peptide_mass_tolerance_upper ${peptide_mass_tolerance_upper} \
        -peptide_mass_tolerance_lower ${peptide_mass_tolerance_lower} \
        -peptide_mass_units ${peptide_mass_units} \
        -isotope_error ${isotope_error} \
        -fragment_bin_tol ${fragment_bin_tol} \
        -fragment_bin_offset ${fragment_bin_offset} \
        -theoretical_fragment_ions ${theoretical_fragment_ions} \
        -static_modifications "${static_modifications}"
    """
}

/*
 * Generates a decoy database from the given FASTA file
 *
 * @return FASTA with decoys
 */
process generate_decoy_database {
    label 'mcquac_image'

    cpus 2
    memory "8.GB"

    publishDir path: { "${fasta_output_folder}" }, mode: 'copy', pattern: "*.fasta", enabled: "${store_decoy_fasta}"  // Publish the FASTA file, which was used for the search

    input:
    path fasta
    val decoy_method
    val store_decoy_fasta
    val fasta_output_folder

    output:
    path "${fasta.baseName}_with_decoys.fasta"

    script:
    """
    DecoyDatabase -in ${fasta} -out ${fasta.baseName}_with_decoys.fasta -method ${decoy_method} -decoy_string DECOY_
    """
}


/*
 * Identifies peptides in MS/MS spectra using Comet
 *
 * @return mzid Path to mzid file
 */
process comet_search {
    label 'comet_image'
    
    cpus { comet_cpu }
    memory { comet_mem }
    
    input:
    path mzml
    path input_fasta
    path comet_params_file
    val comet_cpu
    val comet_mem

    output:
    path "${mzml.baseName}.mzid", emit: mzids
    path input_fasta, emit: fasta_file

    // Run Comet search and then fix potential issues in the mzid file when no identifications were found
    script:
    """
    comet -P${comet_params_file} -D${input_fasta} ${mzml}

    if ! grep -q '<SpectrumIdentificationResult' "${mzml.baseName}.mzid"; then
        echo "fixing"
        sed -i 's;</SpectrumIdentificationResult>;<SpectrumIdentificationResult />;g' "${mzml.baseName}.mzid"
    fi

    """
}
