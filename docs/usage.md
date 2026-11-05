# nf-core/macproqc: Usage

## :warning: Please read this documentation on the nf-core website: [https://mpc-bioinformatics.github.io/macproqc/usage/](https://mpc-bioinformatics.github.io/macproqc/usage/)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

This page describes how to prepare your input data and run the nf-core/macproqc pipeline. For general nf-core and Nextflow concepts (profiles, resource configuration, reproducibility) please also refer to the [nf-core documentation](https://nf-co.re/docs/running/run-pipelines).

## Samplesheet input

You will need to create a samplesheet with information about the MS runs you would like to analyse. Use `--input` to specify its location. The file must be a comma-separated CSV with a header row.

```bash
--input '[path to samplesheet file]'
```

### Samplesheet format

```csv title="samplesheet.csv"
sample,spectrum_file
msrunone,/path/to/data/CONTROL_REP1.raw
msruntwo,/path/to/data/CONTROL_REP2.d.tar.gz
msrunthree,/path/to/data/CONTROL_REP3.mzML.gz
```

| Column          | Description                                                                             |
| --------------- | --------------------------------------------------------------------------------------- |
| `sample`        | Custom sample name. Cannot contain spaces.                                              |
| `spectrum_file` | Full path to the mass spectrometry data file. See supported formats in the table below. |

### Supported input formats

The pipeline auto-detects the instrument vendor and compression type from the file extension:

| Extension                              | Vendor        | Notes                                       |
| -------------------------------------- | ------------- | ------------------------------------------- |
| `.raw`                                 | Thermo Fisher | Converted via ThermoRawFileParser           |
| `.d`                                   | Bruker        | Converted via tdf2mzml                      |
| `.mzML`                                | Any           | Used directly, no conversion needed         |
| `.[raw,d,mzML].[zip, tar.gz, tar.bz2]` | See above     | Will be decompressed and processed as above |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Protein sequence database

The pipeline requires a protein sequence database in FASTA format for peptide identification. A single database is applied globally to all MS runs listed in the samplesheet. If the identification is skipped (see [Running only parts of the pipeline](#running-only-parts-of-the-pipeline)), no database is needed.

```bash
--fasta '[path to FASTA file]'
```

Accepted file extensions: `.fasta`, `.fa`, `.fas`, `.faa` — optionally gzip-compressed (`.gz`).

### Decoy database generation

By default, the pipeline appends reversed decoy sequences to the target database using the OpenMS `DecoyDatabase` tool, producing a combined target-decoy FASTA. This is required for FDR estimation during peptide identification.

| Parameter                 | Default   | Description                                                                                                                                                                                  |
| ------------------------- | --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--skip_decoy_generation` | `false`   | Skip decoy generation. Use this when your FASTA already contains decoy sequences (labelled with the string defined by `--decoy_string`).                                                     |
| `--decoy_string`          | `DECOY_`  | String that is combined with the accession of the protein identifier to indicate a decoy protein. Used also later by tools to estimate the FDR.                                              |
| `--decoy_method`          | `reverse` | Method by which decoy sequences are generated from target sequences. Valid parameters are 'reverse' and 'shuffle'.                                                                           |
| `--save_decoy_database`   | `false`   | Write the generated target-decoy FASTA to `<outdir>/decoy_database/`. Useful for archiving the exact database used in a run or for reuse in a subsequent run with `--skip_decoy_generation`. |

> [!TIP]
> If you run the pipeline repeatedly on the same dataset, generate the decoy database once with `--save_decoy_database`, then pass that file via `--fasta` and set `--skip_decoy_generation` in subsequent runs to avoid redundant processing.

## Peptide identification

The pipeline searches spectra against the (target-decoy) database using [Comet](https://uwpr.github.io/Comet/). Search parameters are configured through a combination of pipeline parameters and an optional Comet parameter template file.

### Comet configuration template

Comet is configured via a `.params` file. The pipeline ships a built-in default template (`assets/default_configs/comet.params`) that is used when no custom file is provided. You can supply your own template to control settings not exposed as pipeline parameters (e.g. enzyme, missed cleavages, ion series):

```bash
--comet_config_template '/path/to/your/comet.params'
```

The template is adjusted at runtime: output format flags, mass tolerances, fragment ion settings, and modifications are overwritten with the values of the pipeline parameters listed below. All other settings in the template are preserved as-is.

### Search parameters

The following parameters are overwritten in the Comet config at runtime:

| Parameter                        | Default | Description                                                                                                  |
| -------------------------------- | ------- | ------------------------------------------------------------------------------------------------------------ |
| `--peptide_mass_tolerance_upper` | `5.0`   | Upper bound of precursor mass tolerance (`peptide_mass_tolerance_upper`)                                     |
| `--peptide_mass_tolerance_lower` | `-5.0`  | Lower bound of precursor mass tolerance (`peptide_mass_tolerance_lower`). Usually negative.                  |
| `--peptide_mass_units`           | `2`     | Units for precursor mass tolerance: `0` = amu, `1` = mmu, `2` = ppm                                          |
| `--isotope_error`                | `2`     | Isotope offset range: `0` = off, `1` = 0/1 (C13), `2` = 0/1/2, `3` = 0/1/2/3, `4` = -1/0/1/2/3, `5` = -1/0/1 |
| `--fragment_bin_tol`             | `0.02`  | Fragment ion bin width in amu (`fragment_bin_tol`)                                                           |
| `--fragment_bin_offset`          | `0.0`   | Fragment ion bin offset, 0.0–1.0 (`fragment_bin_offset`)                                                     |
| `--theoretical_fragment_ions`    | `0`     | `0` = use flanking peaks, `1` = monoisotopic peak only                                                       |
| `--psms_per_spectrum`            | `5`     | Number of PSM candidates reported per spectrum (`num_output_lines`)                                          |

The pipeline always enforces the following settings regardless of what is set in the template:

- Output format is set to mzIdentML (all other Comet output formats are disabled).
- Internal decoy search is disabled (`decoy_search = 0`); decoys are handled externally by the pipeline.
- `max_duplicate_proteins = -1` to report all protein matches per PSM.
- `equal_I_and_L = 0` (isoleucine and leucine are treated as distinct).

### Peptide modifications

Variable and static modifications can be set via pipeline parameters. The syntax matches the Comet params file format; multiple entries are separated by a semicolon. The modification string must be enclosed in single quotes when passed on the command line.

| Parameter                  | Default | Description                                                                                                                      |
| -------------------------- | ------- | -------------------------------------------------------------------------------------------------------------------------------- |
| `--variable_modifications` | `""`    | Variable modifications to set in the Comet params file (`variable_modXX = ...` entries)                                          |
| `--static_modifications`   | `""`    | Static modifications to set in the Comet params file (`add_XX_... = ...` entries)                                                |
| `--label_modifications`    | `""`    | Label-specific static modifications for a second, parallel labelled search (see below). Same syntax as `--static_modifications`. |

> [!TIP]
> When passing modification strings on the command line, wrap the entire value in double quotes and the Comet-format string in single quotes:
>
> ```bash
> --variable_modifications "'variable_mod01 = 15.9949 M 0 3 -1 0 0 0.0;variable_mod02 = 79.966331 STY 0 3 -1 0 0 97.976896'"
> --static_modifications "'add_C_cysteine = 57.021464'"
> ```
>
> Only modifications already present in the template can be overwritten this way. To add new modification slots, provide a custom `--comet_config_template`.

### Labelled search (for spike-ins)

When `--label_modifications` is set to a non-empty string, the pipeline runs a **second, parallel Comet search** alongside the standard unlabelled search. The labelled search uses the same parameters as the unlabelled search, but the `--label_modifications` entries are appended to `--static_modifications` to account for stable isotope labels (e.g. heavy lysine and arginine).

This labelled search is only used for the [spike-ins](#spike-ins): its identifications help to find the labelled spike-in peptides for the XIC extraction. It is not used for any other QC metric, and it is not mandatory for the spike-in metrics. The point of a separate search is that the labels can be searched as static modifications only here, instead of adding them as variable modifications to the general identification, which would enlarge the search space and affect the identification based QC metrics.

Example for spike-ins labelled with heavy lysine (K+8) and arginine (R+10):

```bash
--static_modifications "'add_C_cysteine = 57.021464'"
--label_modifications  "'add_K_lysine = 8.014199;add_R_arginine = 10.008269'"
```

Each spectrum file produces two result files: one from the unlabelled search and one from the labelled search. Both are written to the `comet/` output directory.

### Spike-ins

Specific quality control samples may contain spike-in or reference peptides. If a spike-in table is given, the extracted ion chromatograms (XICs) of these peptides are extracted from the Thermo `.raw` and Bruker `.d` files, and their maximum intensity, retention time and number of identifications are added to the QC metrics and the visualization.

```bash
--spike_ins_table '[path to spike-ins CSV file]'
```

The spike-ins are given in a comma-separated table with the following format (an example is provided with the pipeline in `assets/default_configs/spike_ins_whole_proforma.csv`):

```csv title="spike_ins.csv"
name,sequence,mz,RT,mz-tol,mz-tol-unit,rt-tol
SPIKE1,GEPAAAAAPEAGASPVEK[+8.014199]/2,815.9118,-1,10,ppm,300
SPIKE2,NLVVGDETTSSLR[+10.008269]/2,700.8664,1830,10,ppm,300
```

| Column        | Description                                                                                                                                                           |
| ------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `name`        | Unique identifier for the spike-in peptide.                                                                                                                           |
| `sequence`    | ProForma peptidoform sequence of the spike-in, including any modifications and the charge (e.g. `/2`).                                                                |
| `mz`          | Expected mass-to-charge ratio (m/z) of the spike-in.                                                                                                                  |
| `RT`          | Expected retention time of the spike-in in seconds, or `-1` if it is not known.                                                                                       |
| `mz-tol`      | Mass tolerance around `mz` used to extract the XIC (a number, e.g. `10`).                                                                                             |
| `mz-tol-unit` | Unit of `mz-tol`: `ppm` (recommended). Thermo `.raw` files also support `amu`, `mmu` and `da`, while the XIC extraction from Bruker `.d` folders only supports `ppm`. |
| `rt-tol`      | Retention time tolerance in seconds around the (expected or identified) retention time, used to extract the XIC.                                                      |

The table is validated at the start of the pipeline (whenever it is used by the metrics extraction or the visualization): all columns are required, `name` must be unique, `mz` and `mz-tol` must be greater than 0, `RT` must be `-1` or a value of at least 0, `rt-tol` must be at least 0 and `mz-tol-unit` must be one of `ppm`, `amu`, `mmu` or `da`.

The retention time window used to extract the XIC of each spike-in is determined as follows:

1. If the spike-in was identified in the run, the window is centred on the mean retention time of its identifications. The half window size is `rt-tol`, or half of the retention time range of the identifications, if that is larger.
2. Otherwise, if an `RT` is given, the window is `RT` ± `rt-tol`.
3. Otherwise (`RT` is `-1` and no identification), the XIC is extracted over the whole run.

For this tuning of the retention time windows, the PSMs of the unlabelled **and** (if `--label_modifications` is set) the labelled search are used. The spike-ins are matched to the PSMs by their plain amino acid sequence, so each search contributes the spike-ins it can identify: e.g. unlabelled spike-ins are found in the unlabelled search, while labelled spike-ins are found in the labelled search. The filtering of these PSMs is intentionally more relaxed than the general FDR filtering: all PSMs with a Comet expectation value below 0.01 are used. These identifications are also used to report the number of identifications of each spike-in.

Without identifications (e.g. with `--skip_identification`), only the retention times given in the table are used, and no spike-in is counted as identified.

## FDR filtering and protein inference

The PSM results are processed with [PIA](https://github.com/medbioinf/pia) (PIA - Protein Inference Algorithms) for FDR filtering and peptide/protein inference.

| Parameter                   | Default | Description                                                                                                                                                                    |
| --------------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `--pia_fdr_threshold`       | `0.01`  | FDR threshold used by PIA for the final PSM-, peptide- and protein-level filtering.                                                                                            |
| `--pia_prefilter_threshold` | `0.05`  | FDR threshold for an optional PSM-level pre-filtering pass, run before the final `--pia_fdr_threshold` analysis on the unlabelled search. Set to `0` to disable pre-filtering. |

When `--pia_prefilter_threshold` is greater than `0`, the unlabelled search's PSMs are first filtered on this looser threshold and recompiled, before the final analysis (at `--pia_fdr_threshold`) runs on the reduced PSM set - this speeds up inference on large searches.

If a run has no (FDR-valid) identifications, e.g. because the search finds no PSM at all or no PSM passes the pre-filtering or the final FDR filtering, it is still reported, with 0 PSMs, peptides and protein groups.

Labelled search results (produced when `--label_modifications` is set, mainly used for spike-in peptides) only ever receive a PSM-level export - protein and peptide inference is not meaningful for the labelled search channel, so both are skipped, and no pre-filtering is applied to labelled results either.

## Visualization

The visualization creates plots and a summary table over all MS runs (see the [output documentation](output.md#visualization)). Their formats can be set with:

| Parameter             | Default       | Description                                                                                                      |
| --------------------- | ------------- | ---------------------------------------------------------------------------------------------------------------- |
| `--figure_format`     | `html,plotly` | Format of the figures: `plotly` (plotly JSON files), `html` (standalone interactive HTML files) or `html,plotly` |
| `--output_table_type` | `tsv`         | Format of the summary and loadings tables: `tsv`, `csv` or `xlsx`                                                |

## Running only parts of the pipeline

The pipeline consists of three main steps, which are all executed by default and can each be disabled:

| Parameter                   | Step                                                                                                                                     |
| --------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| `--skip_metrics_extraction` | Metrics extraction: preparation of the spectra, identification and collection of all metrics into one HDF5 file per run (`metrics/hdf5`) |
| `--skip_mzqc_export`        | Export of the metrics into mzQC files (`metrics/mzqc`)                                                                                   |
| `--skip_visualization`      | Visualization of the metrics (`qcvisualization`)                                                                                         |

Additionally, the identification can be skipped within the metrics extraction with `--skip_identification`. Then no identification based metrics are collected, no `--fasta` is needed, and spike-in metrics (if `--spike_ins_table` is given) only use the retention time windows of the table.

If the metrics extraction is skipped, the HDF5 files of an earlier run have to be given with `--input_hdf5_dir` instead of a samplesheet with `--input` (both cannot be combined). This allows e.g. to re-create the visualization with other settings, without repeating the identification:

```bash
nextflow run nf-core/macproqc \
   --skip_metrics_extraction \
   --skip_mzqc_export \
   --input_hdf5_dir ./results/metrics/hdf5 \
   --figure_format html \
   --outdir ./results_html \
   -profile docker
```

All files ending with `.hdf5` in this directory are used, the run names are taken from the file names.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/macproqc \
   --input ./samplesheet.csv \
   --fasta ./database.fasta \
   --outdir ./results \
   -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work/               # Directory containing the Nextflow working files
<OUTDIR>/           # Finished results in the specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other Nextflow hidden files, e.g. history of pipeline runs and old logs
```

### Params file

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify them in a YAML or JSON params file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

```bash
nextflow run nf-core/macproqc -profile docker -params-file params.yaml
```

with `params.yaml`:

```yaml title="params.yaml"
input: "./samplesheet.csv"
fasta: "./database.fasta"
outdir: "./results/"
```

You can also generate parameter files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available — even if the pipeline has been updated since. To make sure that you're running the latest pipeline version, regularly update the cached version:

```bash
nextflow pull nf-core/macproqc
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline.

First, go to the [nf-core/macproqc releases page](https://github.com/mpc-bioinformatics/macproqc/releases) and find the latest pipeline version — numeric only (e.g. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) — e.g. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Apptainer, Conda) — see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` — the order of arguments is important! They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
