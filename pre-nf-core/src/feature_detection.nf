#!/usr/bin/env nextflow

// 
// Workflows for the feature detection
// 

/**
 * Extracts features from peak picked MS1 spectra
 */
workflow get_feature_metrics {
    take:
        mzmls  // MS1 should be peak picked for feature-finding
        mztabfiles
        precursor_tolerance
        precursor_tolerance_unit
        min_charge
        max_charge
        openms_threads
        openms_memory
    
    main:
        // filter out empty spectra and chromatograms
        filtered_mzml = filter_mzml(mzmls, openms_threads, openms_memory)

        // Get features with OpenMS feature finder
        feature_xml = run_feature_finder(
            filtered_mzml,
            precursor_tolerance,
            precursor_tolerance_unit,
            min_charge,
            max_charge,
            openms_threads,
            openms_memory,
        )

        // map identification to features (using mzTab and featureXML)
        runbase_to_featurexml_and_mztab = 
            feature_xml.map {
                file -> tuple(file.name.take(file.name.lastIndexOf('.filtered.featureXML')), file)
            }.join(
                mztabfiles.map{
                    file -> tuple(file.name.take(file.name.lastIndexOf('-piaExport-PSM.mzTab')), file)
                },
                by: 0
            )
        feature_xml_identified = map_features_to_idents(runbase_to_featurexml_and_mztab, openms_threads, openms_memory)

        // Retrieve the actual data and report a HDF5 file
        feature_metrics = get_metrics_from_featurexml(feature_xml_identified, max_charge)
        
    emit:
         feature_metrics
}

process run_feature_finder {
    label 'mcquac_image'

    cpus { openms_threads }
    memory { openms_memory }

    input:
    path mzml
    val precursor_tolerance
    val precursor_tolerance_unit
    val min_charge
    val max_charge
    val openms_threads
    val openms_memory

    output:
    path "${mzml.baseName}.featureXML"

    script:
    """
    tol=${precursor_tolerance}
    tol_unit=${precursor_tolerance_unit}

    if [[ "\${tol_unit}" == "0" ]]; then
        tol_unit="DA"
    elif [[ "\${tol_unit}" == "1" ]]; then
        tol=\$(echo "\${tol} / 1000.0" | bc -l)
        tol_unit="DA"
    elif [[ "\${tol_unit}" == "2" ]]; then
        tol_unit="ppm"
    fi

    FeatureFinderMultiplex -in ${mzml} -out ${mzml.baseName}.featureXML -threads ${openms_threads} \
        -algorithm:labels "" \
        -algorithm:charge "${min_charge}:${max_charge}" \
        -algorithm:spectrum_type centroid \
        -algorithm:mz_tolerance \${tol} \
        -algorithm:mz_unit \${tol_unit}
    """
}

process map_features_to_idents {
    label 'mcquac_image'

    cpus { openms_threads }
    memory { openms_memory }

    input:
    tuple val(runBaseName), path(featurexml), path(ident)
    val openms_threads
    val openms_memory

    output:
    path "${featurexml.baseName}_with_idents.featureXML"

    script:
    """
    convert_mztab_to_idxml.py -mztab ${ident} -out_idxml ${ident.baseName}.idXML

    IDMapper -in ${featurexml} -out ${featurexml.baseName}_with_idents.featureXML -threads ${openms_threads} \
        -id ${ident.baseName}.idXML 
    """
}

process get_metrics_from_featurexml {
    label 'mcquac_image'

    cpus 1
    memory "4.GB"

    input:
    path featurexml
    val max_charge

    output:
    path "${featurexml.name.take(featurexml.name.lastIndexOf('.filtered_with_idents.featureXML'))}-features.hdf5"

    script:
    """
    extract_from_featurexml.py -featurexml ${featurexml} \
        -report_up_to_charge ${max_charge} \
        -out_hdf5 ${featurexml.name.take(featurexml.name.lastIndexOf('.filtered_with_idents.featureXML'))}-features.hdf5 
    """
}

/**
 * Removes chromatrograms and empty peak lists from mzML
 * Necessary to prevent memory issues with FeatureFinder
 */
process filter_mzml {
    label 'mcquac_image'

    cpus { openms_threads }
    memory { openms_memory }

    input:
    path mzml
    val openms_threads
    val openms_memory

    output:
    path "${mzml.baseName}.filtered.mzML"

    script:
    """
    # Filter mzML for MS1 spectra
    FileFilter -in ${mzml} -out ${mzml.baseName}.filtered.mzML -threads ${openms_threads} \
        -peak_options:remove_chromatograms \
        -peak_options:remove_empty \
        -peak_options:sort_peaks \
        -peak_options:zlib_compression true
    """
}