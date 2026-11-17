#!/usr/bin/env nextflow

// Mass centric Quality Control Workflow (McQuaC)
// 
// This is a QC workflow which generates various metrics from Thermo "RAW-files", Bruker "d.-paths" as well as mzML files. 
// Please Read the documentation for more information.


// import nf-schema parameters validation
include { validateParameters } from 'plugin/nf-schema'

// Include all the needed workflows from the sub-workflows
include {convert_raws_to_mzml} from './src/io/raw_file_conversion.nf'
include {identification_with_comet; identification_with_comet as identification_labelled_with_comet} from './src/identification/comet.nf'
include {pia_analysis_full; pia_analysis_psm_only; pia_extract_metrics} from './src/pia.nf'
include {retrieve_spike_ins_information} from './src/retrieve_spike_ins.nf'
include {get_feature_metrics} from './src/feature_detection.nf'
include {get_headers; get_mzml_infos} from './src/metrics/ms_run_metrics.nf'
include {combine_metric_hdf5} from './src/io/combine_metric_hdf5.nf'
include {hdf5_to_mzqc} from './src/mzqc.nf'
include {output_processing_success} from './src/io/output_processing_success.nf'
include {visualization} from './src/visualization.nf'

// MAIN WORKFLOW
workflow {
	// validate the parameters to the schema - if they don't validate, workflow will be stopped, if not in the schema, a warning is thrown
	validateParameters()
			
	if (params.visualize_only) {
		// only visualization
		list_hdf5_files = file(params.main_input_folder + '/*.hdf5')
		visualization(
			list_hdf5_files,
			params.main_outdir,
			params.rt_unit,
			params.output_column_order,
			params.spikein_columns,
			params.output_table_type,
			params.figure_format,
			params.search_spike_ins,
			params.height_barplots,
			params.width_barplots,
			params.height_pca,
			params.width_pca,
			params.height_ionmaps,
			params.width_ionmaps,
			params.visualization_mem, 
			params.spike_ins_table
		)
	} else {
		// complete workflow

		// Retrieve input files
		fasta_file = Channel.fromPath(params.main_fasta_file).first()
		thermo_raw_files = Channel.fromPath(params.main_input_folder + "/*.raw")
		bruker_raw_folders = Channel.fromPath(params.main_input_folder + "/*.d", type: 'dir')
		input_mzml_files = Channel.fromPath(params.main_input_folder + "/*.mzML")
		
		raw_files = thermo_raw_files.concat(bruker_raw_folders)

		// conversion into mzML files
		converted_mzmls = convert_raws_to_mzml(
			thermo_raw_files, 
			bruker_raw_folders, 
			params.file_conversion__thermo_raw_conversion_mem, 
			params.file_conversion__bruker_raw_conversion_cpu, 
			params.file_conversion__bruker_raw_conversion_mem
		)
		mzmls = converted_mzmls.concat(input_mzml_files)
	
		// retrieve mzML Metrics
		mzml_metrics = get_mzml_infos(
			mzmls,
			params.ms_run_metrics__mzml_mem,
			params.base_peak_tic_up_to,
			params.filter_threshold,
			params.report_up_to_charge
		)

		// Identify spectra using Comet
		comet_ids = identification_with_comet(
			mzmls,
			fasta_file,
			params.identification__generate_decoys,
			params.identification__decoy_method,
			params.main_outdir,
			params.identification__store_decoy_fasta,
			params.identification__comet_threads,
			params.identification__comet_mem,
			params.identification__peptide_mass_tolerance_upper,
			params.identification__peptide_mass_tolerance_lower,
			params.identification__peptide_mass_units,
			params.identification__isotope_error,
			params.identification__fragment_bin_tol,
			params.identification__fragment_bin_offset,
			params.identification__theoretical_fragment_ions,
			""	// no labelling searched here
		)

		// Execute protein inference and filter by FDR
		pia_report_files = pia_analysis_full(
			comet_ids.mzids,
			params.identification__pia_threads,
			params.identification__pia_gb_ram,
			params.identification__pia_fdr_threshold,
			params.identification__pia_prefilter
		)
		pia_report_psm_mztabs = pia_report_files
				.toList()
				.transpose()
				.first()
				.flatten()
		pia_extract_csv = pia_extract_metrics(pia_report_files, params.identification__peptides_table, params.main_outdir)

		// search additionally for labelled PSMs
		if (params.search_labelled_spikeins) {
			comet_labelled_ids = identification_labelled_with_comet(
				mzmls,
				fasta_file,
				params.identification__generate_decoys,
				params.identification__decoy_method,
				params.main_outdir,
				params.identification__store_decoy_fasta,
				params.identification__comet_threads,
				params.identification__comet_mem,
				params.identification__peptide_mass_tolerance_upper,
				params.identification__peptide_mass_tolerance_lower,
				params.identification__peptide_mass_units,
				params.identification__isotope_error,
				params.identification__fragment_bin_tol,
				params.identification__fragment_bin_offset,
				params.identification__theoretical_fragment_ions,
				params.identification__label_modifications
			)

			// set the filter to true to count only FDR filtered labelled PSMs - but the FDR is skewed anyways, as the labelling is set to "static"!
			labelled_pia_report_files = pia_analysis_psm_only(
				comet_labelled_ids.mzids,
				false,
				params.identification__pia_threads,
				params.identification__pia_gb_ram,
				params.identification__pia_fdr_threshold
			)
		}

		// extract spike-ins information
		if (params.search_spike_ins) {
			spike_ins_table = Channel.fromPath(params.spike_ins_table).first()

			if (params.search_labelled_spikeins) {
				psm_results = labelled_pia_report_files
			} else {
				psm_results = pia_report_psm_mztabs
			}

			spike_in_metrics = retrieve_spike_ins_information(
				raw_files,
				psm_results,
				spike_ins_table,
				params.max_parallel_xic_extractors_factor
			)
		}
	 
		// Run Feature Finding
		feature_metrics = get_feature_metrics(
			mzmls,
			pia_report_psm_mztabs, 
			params.identification__peptide_mass_tolerance_upper,
			params.identification__peptide_mass_units,
			params.feature_detection__min_charge,
			params.feature_detection__max_charge,
			params.feature_detection__openms_threads,
			params.feature_detection__openms_memory
		)

		// Get Thermo/Bruker specific information from raw_spectra
		custom_header_infos = get_headers(
			thermo_raw_files,
			params.ms_run_metrics__thermo_raw_mem,
			params.ms_run_metrics__thermo_headers,
			bruker_raw_folders,
			params.ms_run_metrics__bruker_raw_mem,
			params.ms_run_metrics__bruker_headers,
			params.ms_run_metrics__bruker_calibrants
		)

		// Concatenate to one merged metric CSV
		hdf5s_per_run = mzml_metrics.map{file -> tuple(file.name.take(file.name.lastIndexOf('-mzml_info.hdf5')), file)}
		if (params.search_spike_ins) {
			hdf5s_per_run = hdf5s_per_run.concat(
				spike_in_metrics.map{file -> tuple(file.name.take(file.name.lastIndexOf('-spikeins.hdf5')), file)}
			)
		}
		hdf5s_per_run = hdf5s_per_run
			.concat(feature_metrics.map{file -> tuple(file.name.take(file.name.lastIndexOf('-features.hdf5')), file)})
			.concat(pia_extract_csv.map{file -> tuple(file.name.take(file.name.lastIndexOf('-pia_extraction.hdf5')), file)})
			.concat(custom_header_infos.map{file -> tuple(file.name.take(file.name.lastIndexOf('-custom_headers.hdf5')), file)})
			.groupTuple()

		combined_metrics = combine_metric_hdf5(hdf5s_per_run, params.main_outdir)

		hdf5_to_mzqc(combined_metrics, params.main_outdir)

		// Visualize the results (and move them to the results folder)
		visualization(
			combined_metrics,
			params.main_outdir,
			params.rt_unit,
			params.output_column_order,
			params.spikein_columns,
			params.output_table_type,
			params.figure_format,
			params.search_spike_ins,
			params.height_barplots,
			params.width_barplots,
			params.height_pca,
			params.width_pca,
			params.height_ionmaps,
			params.width_ionmaps, 
			params.visualization_mem, 
			params.spike_ins_table
		)

		output_processing_success(
			raw_files.concat(input_mzml_files),
			hdf5s_per_run.toList().transpose().first().flatten(),
			params.main_outdir
		)
	}
}
