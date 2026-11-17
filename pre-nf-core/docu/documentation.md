# Documentation for the Nextflow QC workflow

This workflow does quality control for DDA mass spectrometry proteomics data. 
The workflow can be used in two different modes: for ISA runs (is_isa = TRUE) or normal runs (is_isa = FALSE).

Some metrics are based on the following paper:
Bittremieux W, Meysman P, Martens L, Valkenborg D, Laukens K. Unsupervised Quality Assessment of Mass Spectrometry Proteomics Experiments by Multivariate Quality Control Metrics. J Proteome Res. 2016 Apr 1;15(4):1300-7. doi: 10.1021/acs.jproteome.6b00028. Epub 2016 Mar 18. PMID: 26974716.


## Input: 

- raw files (TODO: from where?)
- FASTA file(s) (TODO: from where?)
- Config files for Comet, MSGF+ and PIA (TODO: from where?, are they fixed or dynamically generated from user input?)
- List with spike-ins (only for ISA runs)


## Output:
- A table with quality control metrics is saved to the database.
- Figures for parts of these metrics are created.

## General workflow:
- **Data input**: TODO
- **Raw data metrics**: Raw files are converted to mzML, then different metrics are collected and calculated (e.g., number of MS1 and MS2 spectra, TICs, precuror charges).
- **Peptide identification**: The MS2 spectra in the raw files are identified using the Comet and the MSGF+ search engines using the provided FASTA files.
- **FDR filtering and protein inference**: The results from Comet and MSGF+ are further analysed by PIA (aggregating the results of the two search engines, filter for FDR and performing protein inference). PIA creates a mzTab file, that is further processed.
- **Spike-ins**: Only for ISA: Identifications and XIC for a specific list of spike-ins are searched in the mzTab and the corresponding 
- **Feature finding**: Using OpenMS FeatureFinder to find features and get metrics like total and identified number of features.
- **Saving metrics to database**: The output table is saved to the database.
- **Plotting**: In a separate step, information is retrieved from the data base and plots are created (barplots for various metrics, PCA plots and TIC overlay).

<img src="Workflow_Foto.jpg"  width="600" height="600">

*TODO: make the workflow graph nicer*.

## Details

### 1) **Data input**

### 2) **Raw data metrics**
Sub-workflow: get_mzml_chromatogram_and_more.nf

The raw files are first converted to mzML using the ThermoRawFileParser with peak picking. The Python script extract_data_from_mzml.py uses the pyopenms package to calculate various metrics like number of MS1 and MS2 spectra, information on retention time and TIC, precursor charge states etc.
For details on the different metrics see the QC_columns_overview.csv table.

### 3) **Peptide identification**
Sub-workflows: identification_via_comet.nf and identification_via_msgfplus.nf  

The MS2 spectra in the raw files are identified via the Comet and MSGF+ search engines using the provided FASTA files.

For Comet, a configuration file (example_configurations/high-high.comet.params) is needed, which contains the necessary parameters. 
By default, peptide mass tolerance is set to 5ppm, fragment tolerance as 20ppm, maximum 2 missed cleavages and Trypsin as the enzyme.
Oxidation(M) is set as a variable modification and Cabamidomethylation (C) as a fixed modification.

For MSGF+ (under construction), the configuration file example_configurations/msgfplus_config.txt is used. Raw files have to be converted to MGF files via the ThermoRawFileParser before identification.
By default, instrument type is set to high-resolution, fragment tolerance as 20ppm, unlimited number of missed cleavages and Trypsin as the enzyme.
Oxidation(M) is set as a variable modification and Cabamidomethylation (C) as a fixed modification.

### 4) **FDR filtering and protein inference**
Sub-Workflow: pia.nf

For FDR estimation and protein inference, the CLI version of PIA is used.
The idXML files from the identification steps are used and processed using PIA compilation and analysis.

The PIA results are then given as a piaExport-PSMs.mzTab, piaExport-peptides.csv and piaExport-proteins.mzid.

From these files, also metrics like the number of identified proteins, proteingroups, PSM charges and missed cleavages are calculated (TODO!).

*Are the results of Comet and MSGF+ also combined here?*

### 5) **Spike-ins**
Sub-Workflow: retrieve_spike_ins_thermorawfileparser.nf

For ISA files only, metrics for a set of predefined spike-in peptides can be calculated (see QC_columns_overview.csv).
For this, the spike-in peptides have to be included in the corresponding FASTA file. 
The list of currently used spike-ins is given in spike_ins.csv.
Using ThermoRawFileParser is used to retrieve the XICs within 10ppm m/z window and +/- 3 min retention time window of the expected spike-in values.
From the output of ThermoRawFileParser we aquire the following information for each spike-in:
- FIXED_MPCSPIKEX_PEP_XXX_RT_XXX -> Sum of all peaks measured inside the expected spike-in window (10 ppm m/z-tolerance and +/- 3 min retention time window).
- IDENT_MPCSPIKEX_COUNT -> spectral count for the given peptide sequence.
- IDENT_MPCSPIKEX_DELTA_RT -> difference between expected and measured retention time for the specific spike-in (if peptide was identified).
- IDENT_MPCSPIKEX_PEP_XXX_RT_DELTA -> Sum of XIC for the identified PSMs.

### 6) **Feature Finding**
Sub-workflow. get_features_in_raws.nf

The raw files are converted to mzML (*Why don't we use the mzMLs from step 2?*) using ThermoRawFileParser. Then, the OpenMS FeatureFinderCentroided is used to find the features. The OpenMS IDmapper is used to match the found features with the identifications (mzTabs acquired from PIA). The resulting feature.xml file is then used to calculate further metrics like total number of features, charge state of features etc. (see QC_columns_overview.csv) .


### 7) **Saving results to database**
*TODO!*

### 8) **Plotting**
The visualization is done in Python with the plotly package for interactive plots.

The data are directly read from the csv file provided by the previous workflow step. 
The following parameters can be set:
group (default:None): List of experimental group (comma-separated), will influence the colouring in the PCA plots.
tic_overlay_offset (default: 0): allows offset for the TIC plots (so there is a distance between them and you can see them better).


The data stem from the different tables in the database. 

List of plots/tables for display:

- Table 0: Overview over QC measures. This table shows an overview over the different QC measures (one row per raw file).  
- Figure 1: Number of MS1 and MS2 spectra. This figure shows the total number of MS1 and MS2 spectra measured in each sample.  
- Figure 2: Number of PSMs, peptides and protein groups. This figure shows a barplot of the number of peptide-spectrum matches (PSMs) and peptides, each filtered by a false discovery rate of 1%. Additionally, the protein groups formed by PIA as well as the total number of protein accessions in these groups (ungrouped proteins) are shown.
- Figure 3: Number of features and identified features. This figure shows a barplot of the total number of features and those that were identified in the raw data.
- Figure 4: TIC overlay. This figure shows the overlay of the total ion chromatograms (TICs) of all raw files. On the x-axis the retention time and on the y-axis the intensity of the TIC is given.
- Figure 5: TIC quartiles. This barplot shows the TIC quartiles, so how many % of the retention time was needed to measure the first, second, third or fourth 25% of the total TIC. The similarity of TICs between the different samples can be assessed by this figure.
- Figure 6: MS1 TIC quartiles. This barplot shows the MS1 TIC quartiles, so how many % of the retention time was needed to measure the first, second, third or fourth 25% of the TIC restricted to MS1. The similarity of TICs between the different samples can be assessed by this figure.
- Figure 7: MS2 TIC quartiles. This barplot shows the MS2 TIC quartiles, so how many % of the retention time was needed to measure the first, second, third or fourth 25% of the TIC restricted to MS2. The similarity of TICs between the different samples can be assessed by this figure.
- Figure 8: Precursor charge states. This barplot show the relative distribution of charge states (1-5 or more and unknown) of the precursor ions, so before any identification. 
- Figure 9: PSM charge states. This barplot show the relative distribution of charge states (1-5 or more and unknown) of the ions with a peptide-spectrum match, so after identification. 
- Figure 10: PSM missed cleavages. This barplot show the relative distribution of missed cleavages (0-3 or more) of the identified peptides. 
- Figure 11: PCA on all data (a) + loading plot (b). In the PCA each dot stands for one raw file. For calculation of the PCA a set of QC measures before and after identification is used. In the loading plots each dot stands for one of these Qc measures. Important QC measures for the shape of the PCA plot can be derived by looking at the upper, lower and left and right borders. 
- Figure 12: PCA on raw data (a) + loading plot (b). In the PCA each dot stands for one raw file. For calculation of the PCA a set of QC measures is used, which can be derived directly from the raw data (before any identification). In the loading plots each dot stands for one of these Qc measures. Important QC measures for the shape of the PCA plot can be derived by looking at the upper, lower and left and right borders.
- Figure 13: Ion Maps. For each raw file an ion map is generated. The x-axis shows the retention time and the y-axis the m/z. The colouring show the MS2 TIC.
- Figure 14: Pump pressure. This line plot show the pump pressure (y-axis) depending on the retention time (x-axis). 
- Figure 15: Various additional plots. These line plots show different further variables depending on the retention time, for example ion injection time or lock mass correction.



