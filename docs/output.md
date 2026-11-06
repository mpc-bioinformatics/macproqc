# nf-core/macproqc: Output

## Introduction

This document describes the output produced by the pipeline. All paths are relative to the top-level results directory set with `--outdir`.

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes mass spectrometry data through the following steps:

- [Spectra preparation](#spectra-preparation) - decompression and conversion to mzML
- [Peptide identification](#peptide-identification) - decoy database generation and database search with Comet
- [FDR filtering and protein inference](#fdr-filtering-and-protein-inference) - FDR filtering and PIA protein groups
- [QC metrics](#qc-metrics) - all QC metrics of an MS run, collected into one HDF5 file
- Feature finding - isotope feature detection and feature based QC metrics _(migrating, not yet part of the pipeline)_
- [mzQC output](#mzqc-output) - standardised QC metric export
- [Visualization](#visualization) - graphical visualizations of QC metrics using interactive plotly plots
- [Pipeline information](#pipeline-information) - reports from Nextflow itself

Depending on the skipped steps (see the [usage documentation](usage.md#running-only-parts-of-the-pipeline)), some of these outputs are not created. When the metrics extraction is skipped and the HDF5 files of an earlier run are given with `--input_hdf5_dir`, only the mzQC and/or visualization outputs are written.

---

## Spectra preparation

Raw vendor files are decompressed (if necessary) and converted to the open mzML format. Thermo `.raw` files are converted using ThermoRawFileParser; Bruker `.d` folders are converted using tdf2mzml. Files already in mzML format are passed through after optional decompression.

<details markdown="1">
<summary>Output files</summary>

- `mzmls/`
  - `<sample>.mzML` - mzML files converted from Thermo `.raw` and Bruker `.d` input files (input files already given as mzML are not copied here)

</details>

---

## Peptide identification

Peptide spectrum matches (PSMs) are produced by searching the converted mzML files against a protein sequence database using the [Comet](https://uwpr.github.io/Comet/) search engine.

### Decoy database

Before the database search, the pipeline generates a target-decoy FASTA by appending decoy sequences (reversed by default) to the input database using the OpenMS `DecoyDatabase` tool. This step can be skipped with `--skip_decoy_generation` if the supplied FASTA already contains decoys.

The decoy database file is an intermediate used downstream for FDR estimation. It is only written to the output directory when `--save_decoy_database true` is set.

<details markdown="1">
<summary>Output files</summary>

- `decoy_database/`
  - `*_decoy.fasta` - combined target-decoy protein sequence database. Only present when `--save_decoy_database true` is set.

</details>

### Comet search configuration

Before the database search, the pipeline generates a Comet parameter file by taking a template (either the built-in default at `assets/default_configs/comet.params` or a user-supplied file via `--comet_config_template`) and overwriting the search settings (mass tolerances, fragment ion parameters, modifications, and output format flags) with the pipeline parameter values.

The adjusted parameter file is written next to the search results (see below), to document the settings of each search.

### Comet PSMs

Each mzML file is searched against the target-decoy database using Comet and the adjusted parameter file described above. Results are written as mzIdentML (`.mzid`) files.

When `--label_modifications` is set, two searches are run per spectrum file - one unlabelled and one labelled - producing a pair of mzid files per input. The filename encodes the search type:

- `<sample>-unlabelled.mzid` - standard search using `--static_modifications`
- `<sample>-labelled.mzid` - search with `--label_modifications` merged into the static modifications

<details markdown="1">
<summary>Output files</summary>

- `comet/`
  - `*-unlabelled.mzid` - PSM results from the unlabelled search
  - `*-unlabelled.comet.params` - the adjusted Comet parameter file used for the unlabelled search
  - `*-labelled.mzid` - PSM results from the labelled search (only present when `--label_modifications` is set)
  - `*-labelled.comet.params` - the adjusted Comet parameter file used for the labelled search (only present when `--label_modifications` is set)

</details>

---

## FDR filtering and protein inference

The PSM results are compiled into a [PIA](https://github.com/medbioinf/pia) (PIA - Protein Inference Algorithms) intermediate file, filtered at `--pia_fdr_threshold` FDR, and used for peptide and protein inference. When `--pia_prefilter_threshold` is greater than `0` (the default is `0.05`), an additional pre-filtering pass first reduces the unlabelled PSM set on this looser FDR threshold before the final analysis is run, which speeds up inference on large searches. Labelled search results (produced when `--label_modifications` is set, mainly used for labelled spike-ins) only ever receive a PSM-level export - protein and peptide inference is not meaningful for the labelled search channel, so both are skipped, and no pre-filtering is applied to labelled results either.

<details markdown="1">
<summary>Output files</summary>

- `pia/`
  - `*-unlabelled.piaExport-PSMs.mzTab` - PSM-level results of the unlabelled search, filtered at `--pia_fdr_threshold`
  - `*-unlabelled.piaExport-peptides.csv` - inferred peptides for the unlabelled search
  - `*-unlabelled.piaExport-proteins.mzTab` - inferred protein groups for the unlabelled search
  - `*-labelled.piaExport-PSMs.mzTab` - PSM-level results of the labelled search (only present when `--label_modifications` is set); no peptide- or protein-level export is produced for the labelled search

</details>

---

## QC metrics

For each MS run, the QC metrics are extracted by several steps and then combined into one HDF5 file. The metrics are annotated with their (PSI-MS) controlled vocabulary accessions, where available.

- **mzML metrics** (from the mzML file, using pyOpenMS): numbers of MS1 and MS2 spectra, retention time range, MS1/MS2 total ion current (TIC) over time, TIC quartiles, spectra density quantiles, base peak and TIC maxima, precursor charge fractions, acquisition frequencies and the data for the MS1 maps.
- **Instrument metrics** (from the vendor file): headers and logs recorded by the instrument, e.g. pump pressure, ion injection times, temperatures and lock mass corrections for Thermo `.raw` files, or frame headers and calibrants for Bruker `.d` folders. These are not available for runs given as mzML.
- **Identification metrics** (from the unlabelled PIA results): numbers of PSMs, peptides, protein groups and accessions, PSM charge fractions, missed cleavage fractions and the precursor mass error (in ppm) of the PSMs. Not available with `--skip_identification`.
- **Spike-in metrics** (only with `--spike_ins_table`, for Thermo and Bruker runs): for each spike-in/reference peptide of the table, the maximum intensity and retention time of its extracted ion chromatogram (XIC) and its number of identifications (from the labelled and unlabelled search, if available).

The per-step HDF5 files are intermediates, only the combined file is written to the output directory.

<details markdown="1">
<summary>Output files</summary>

- `metrics/hdf5/`
  - `<sample>.hdf5` - all QC metrics of the MS run. These files can be used as `--input_hdf5_dir` in a later run, e.g. to repeat the visualization with other settings.

</details>

---

## mzQC output

The QC metrics of each MS run are exported into the [mzQC](https://hupo-psi.github.io/mzQC/) format of the HUPO-PSI, a JSON-based standard format for quality control metrics. This step can be skipped with `--skip_mzqc_export`.

<details markdown="1">
<summary>Output files</summary>

- `metrics/mzqc/`
  - `<sample>.mzqc` - the QC metrics of the MS run in mzQC format

</details>

---

## Visualization

The main output of the visualization are interactive plots over all MS runs (as plotly JSON and/or HTML, depending on `--figure_format`). Also, a summary table is generated (as csv, tsv or xlsx, depending on `--output_table_type`). This step can be skipped with `--skip_visualization`.

Some plots are only created, if the respective metrics are available (e.g. instrument headers, Bruker calibrants or identification based metrics). The feature based metrics (number of features and identified features) are not yet migrated to this pipeline, so they are currently not available in the plots.

<details markdown="1">
<summary>Output files</summary>

- `qcvisualization/`
  - `00_table_summary.[csv,tsv,xlsx]` - summary table containing the QC metrics per MS run.
  - `99_hdf5_feature_table.[csv,tsv,xlsx]` - overview of the QC metric keys and names stored in the HDF5 files.
  - `metric_registry.json` - registry of the QC metric types found in the HDF5 files, used to create the plots and tables.
  - `fig01_barplot_MS1_MS2.[plotly.json,html]` - barplot of the number of MS1 and MS2 spectra per MS run.
  - `fig02_barplot_PSMs_peptides_proteins.[plotly.json,html]` - barplot of the number of PSMs, peptides, protein groups and accessions per MS run.
  - `fig03_barplot_features.[plotly.json,html]` - barplot of the number of features and identified features per MS run (feature metrics are not yet migrated).
  - `fig04_MS1_TIC_overlay.[plotly.json,html]` - lineplot of the MS1 total ion chromatogram (TIC) over the retention time for each MS run.
  - `fig05_barplot_TIC_quartiles.[plotly.json,html]` - TIC quartiles are the percentage of retention time that is needed to obtain the first, second, third or fourth 25% of the TIC. These quartiles are visualized as stacked barplots.
  - `fig06_barplot_MS1_TIC_quartiles.[plotly.json,html]` - same as fig05, but only using MS1-level data.
  - `fig07_barplot_MS2_TIC_quartiles.[plotly.json,html]` - same as fig05, but only using MS2-level data.
  - `fig08_barplot_precursor_charge.[plotly.json,html]` - stacked barplot of the precursor charge fractions (also includes unidentified precursors).
  - `fig09_barplot_PSM_charge.[plotly.json,html]` - stacked barplot of the PSM charge fractions (only identified precursors).
  - `fig10_barplot_PSM_missedcleavages.[plotly.json,html]` - stacked barplot of the missed cleavage fractions (only identified).
  - `fig11a_PCA_raw.[plotly.json,html]` - principal component analysis plot (PCA plot). The PCA is calculated on a set of QC metrics that are all directly obtainable from the raw data without a peptide identification step: RT_range, nr_MS1, nr_MS2, accumulated_MS1_TIC, accumulated_MS2_TIC, base_peak_intensity_max, total_ion_current_max, MS2_prec_charge_fraction, RT_MS1_quantiles, RT_MS2_quantiles, RT_TIC_quantiles, MS1_freq_max, MS2_freq_max, MS1_density_quantiles, MS2_density_quantiles, MS1_TIC_change_quantiles, MS1_TIC_quantiles.
  - `fig11b_Loadings_raw.[plotly.json,html]` - loadings plot belonging to the PCA in fig11a. The loadings show which QC metrics have the highest weight in each direction and can help to interpret groups and outliers in the PCA.
  - `fig11c_table_loadings_raw.[csv,tsv,xlsx]` - table of the loadings shown in fig11b.
  - `fig12a_PCA_all.[plotly.json,html]` - PCA plot calculated on the QC metrics used for fig11a plus the metrics obtained during the peptide identification (nr_PSMs, nr_peptides, nr_protein_groups, nr_accessions, PSM_charge_fractions, PSM_missed_cleavage_counts and, once migrated, the feature metrics).
  - `fig12b_Loadings_all.[plotly.json,html]` - loadings plot belonging to the PCA in fig12a.
  - `fig12c_table_loadings_all.[csv,tsv,xlsx]` - table of the loadings shown in fig12b.
  - `fig13_MS1_map/` - folder containing an MS1 map (PNG) for each MS run.
  - `fig14_Pump_pressure.[plotly.json,html]` - lineplot of the pump pressure over the retention time for all MS runs (only if recorded by the instrument).
  - `fig15_PSM_error_boxplots.[plotly.json,html]` - boxplots of the PSM precursor mass error (in ppm) for all MS runs.
  - `fig16_additional_headers/` - folder with lineplots over the retention time for further metrics recorded by the instrument (e.g. temperatures, lock mass correction).
  - `fig17_BRUKER_calibrants/` - folder with lineplots over the retention time for the calibrants of Bruker instruments (only for Bruker runs).

</details>

---

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report_<timestamp>.html` - Nextflow execution report with run statistics and per-process resource usage.
  - `execution_timeline_<timestamp>.html` - timeline of all process executions.
  - `execution_trace_<timestamp>.txt` - tab-separated trace file with per-task resource metrics.
  - `pipeline_dag_<timestamp>.html` - directed acyclic graph (DAG) of the pipeline workflow.
  - `nf_core_macproqc_software_versions.yml` - versions of all software used in the pipeline run.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
