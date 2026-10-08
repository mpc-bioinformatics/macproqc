#  pipeline parameters



## Main parameters

Main parameters for McQuaC

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `main_input_folder` | Folder containing the spectra files to be analysed (Thermo raw, Bruker .d or mzML), or hdf5 files if running only the visualization. All applicable files in this folder will be processed. | `string` | ./ | True |  |
| `main_outdir` | Output-Directory of the result files | `string` | ./results |  |  |
| `main_fasta_file` | The FASTA file used for identification | `string` |  | True |  |

## Visualization parameters

Parameters for visualization script

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `visualize_only` | If true, only the visualization part of teh workflow will be executed. | `boolean` |  |  |  |
| `rt_unit` | Unit of the retention time, either sec for seconds or min for minutes. | `string` | sec |  |  |
| `output_column_order` | Order of columns in the output table | `string` | '' |  |  |
| `spikein_columns` | Columns of the spike-in dataframes that should end up in the result table | `string` | Maximum_Intensity,RT_at_Maximum_Intensity,PSMs,Delta_to_expected_RT |  |  |
| `output_table_type` | Type of the output table, either csv or xlsx | `string` | csv |  |  |
| `height_barplots` | Height of the barplots in pixels | `integer` | 700 |  |  |
| `width_barplots` | Width of the barplots in pixels, 0 = flexible width | `integer` | 0 |  |  |
| `height_pca` | Height of the PCA plot in pixels | `integer` | 1000 |  |  |
| `width_pca` | Width of the PCA plot in pixels | `integer` | 1000 |  |  |
| `height_ionmaps` | Height of the ion maps in inches | `integer` | 10 |  |  |
| `width_ionmaps` | Width of the ion maps in inches | `integer` | 10 |  |  |
| `visualization_mem` | Memory used for the visualization step | `string` | 1.GB |  |  |

## Spike in parameters



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `spike_ins_table` | The information about spike-ins | `string` | ${projectDir}/example_configurations/spike_ins.csv |  |  |
| `search_spike_ins` | Parameter to check if a specific identification and feature extraction for spike in peptides / proteins should be performed. NOTE: FASTA has to contain the spike ins too! | `boolean` | True |  |  |
| `search_labelled_spikeins` | Perform a special identification and look for labelled peptides | `boolean` | True |  |  |

## Metric extraction parameters

Parameters used for the extraction of metrics

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `ms_run_metrics__bruker_raw_mem` | Memory used for the extraction of metrics from Bruker headers | `string` | 1.GB |  |  |
| `ms_run_metrics__thermo_raw_mem` | Memory used for the extraction of metrics from Thermo headers | `string` | 10.GB |  |  |
| `ms_run_metrics__mzml_mem` | Memory used for the extraction of metrics from mzML files | `string` | 10.GB |  |  |
| `ms_run_metrics__bruker_headers` | Set if you want to extract specific headers from Bruker measurements, otherwise the default is used. | `string` |  |  |  |
| `ms_run_metrics__bruker_calibrants` | Set if you want to extract specific calibrants in Bruker raw measurements, otherwise the default is used (622.0290 m/z, 922.009798 m/z and 1221.990637 m/z with a 10 m/z and 0.1 1/k0 tolerance). Have a look into the corresponding python script for the headers. | `string` |  |  |  |
| `ms_run_metrics__thermo_headers` | Set if you want to extract specific headers from Thermo measurements, otherwise the default is used. | `string` |  |  |  |
| `base_peak_tic_up_to` | Retrieve the Basepeak Intensity Max and the Total Ion Current from minute 0 up to the given number in minutes. Defaults to 105 (minutes). | `integer` | 105 |  |  |
| `filter_threshold` | Threshold for the MS1 peaks, to be included in the output file. Defaults to 0.00001 (0.001%) of the highest overall MS1 peak. Values lower will be disregarded. | `number` | 1e-05 |  |  |
| `report_up_to_charge` | Upper limit of range to be reported in a csv table for the charge, defaults to 5. | `integer` | 5 |  |  |
| `max_parallel_xic_extractors_factor` | Factor for the number of maximum forks for the XIC extractor processes, calculated by dividing the number of available processors by this number. In fact, the number of required CPUs per task is set to the given value. | `integer` | 2 |  |  |
| `feature_detection__openms_threads` | Allowed number of threads / CPUs for the OpenMS feature detection | `integer` | 8 |  |  |
| `feature_detection__openms_memory` | Allowed memory for the OpenMS feature detection tasks | `string` | 8.GB |  |  |
| `feature_detection__min_charge` | Minimum charge for the feature detection | `integer` | 2 |  |  |
| `feature_detection__max_charge` | Maximum charge for the feature detection | `integer` | 5 |  |  |

## Conversion parameters



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `file_conversion__bruker_raw_conversion_mem` | Memory for the tdf2mzml to convert .d into mzML | `string` | 5.GB |  |  |
| `file_conversion__bruker_raw_conversion_cpu` | Number of CPUs / threads for the bruker raw file conversion | `integer` | 4 |  |  |
| `file_conversion__thermo_raw_conversion_mem` | Memory for the Thermo Raw File Parser to convert .raw into mzML | `string` | 10.GB |  |  |

## Identification parameters



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `identification__comet_threads` | CPUs / threads used for the identification with Comet | `integer` | 8 |  |  |
| `identification__comet_mem` | Memory used for the identification with Comet | `string` | 10.GB |  |  |
| `identification__generate_decoys` | If true, decoys are generated before starting the spectrum identification. Otherwise, decoys should already be in the FASTA file | `boolean` | True |  |  |
| `identification__store_decoy_fasta` | If true, the generated decoy DB will be stored in the output folder | `boolean` | False |  |  |
| `identification__decoy_method` | Applied method for the generation of decoys, either "reverse" or "shuffle" | `string` | shuffle |  |  |
| `identification__peptide_mass_tolerance_upper` | Comet parameter: peptide_mass_tolerance_upper, upper bound of the precursor mass tolerance | `number` | 5 |  |  |
| `identification__peptide_mass_tolerance_lower` | Comet parameter: peptide_mass_tolerance_lower, lower bound of the precursor mass tolerance; USUALLY NEGATIVE TO BE LOWER THAN 0 | `number` | -5 |  |  |
| `identification__peptide_mass_units` | Comet parameter: peptide_mass_units, 0=amu, 1=mmu, 2=ppm | `integer` | 2 |  |  |
| `identification__isotope_error` | Comet parameter: isotope_error, 0=off, 1=0/1 (C13 error), 2=0/1/2, 3=0/1/2/3, 4=-1/0/1/2/3, 5=-1/0/1 | `integer` | 2 |  |  |
| `identification__fragment_bin_tol` | Comet parameter: fragment_bin_tol, binning to use on fragment ions | `number` | 0.02 |  |  |
| `identification__fragment_bin_offset` | Comet parameter: fragment_bin_offset, offset position to start the binning (0.0 to 1.0) | `number` | 0 |  |  |
| `identification__theoretical_fragment_ions` | Comet parameter: theoretical_fragment_ions, 0=use flanking peaks, 1=M peak only | `integer` | 0 |  |  |
| `identification__label_modifications` | The labels encoded in Comet static labels. Only used during the special identification for labelled spike ins if "search_labelled_spikeins" is true. | `string` | add_K_lysine = 8.014199;add_R_arginine = 10.008269 |  |  |
| `identification__pia_threads` | Allowed number of CPUs / threads for PIA | `integer` | 8 |  |  |
| `identification__pia_gb_ram` | Allowed RAM to be used by PIA, in GB | `integer` | 16 |  |  |
