# nf-core/macproqc

**nf-core/macproqc** (**Ma**ss **C**entric **Pro**teomics **Q**uality **C**ontrol) is a bioinformatics pipeline for comprehensive quality control of mass spectrometry-based proteomics experiments. It accepts raw instrument files from Thermo Fisher (`.raw`) and Bruker (`.d`) instruments as well as pre-converted mzML files and produces a rich set of QC metrics, interactive visualisations, and standardised [mzQC](https://hupo-psi.github.io/mzQC/) output.

The pipeline performs the following main steps:

- **Spectra preparation** - raw vendor files are decompressed and converted to the open mzML format using ThermoRawFileParser (Thermo) or tdf2mzml (Bruker).
- **Peptide identification** - MS2 spectra are searched with Comet against a user-supplied FASTA database (decoys are generated with OpenMS), without and, if label modifications are given, additionally with (static) label information.
- **FDR filtering and protein inference** - PSMs are filtered at a configurable FDR (default 1 %), and peptides and protein groups are inferred with PIA - Protein Inference Algorithms.
- **QC metrics extraction** - per MS run, all QC metrics are collected into one HDF5 file:
  - from the mzML files (pyOpenMS): MS1/MS2 spectrum counts, TIC and TIC quartiles, precursor charge distributions, retention-time coverage, MS1 maps and more
  - from the vendor files: instrument headers and logs, e.g. pump pressure, temperatures and lock mass corrections (Thermo) or calibrants (Bruker)
  - from the identification results: numbers of PSMs, peptides and protein groups, PSM charges, missed cleavages and precursor mass errors
  - optionally for spike-in/reference peptides: XIC-based maximum intensity, retention time and number of identifications
- **QC metrics for peptide features** _(migrating)_ - OpenMS identifies isotope features; found features are mapped to identifications to yield identification rates and QC metrics are extracted.
- **mzQC output** - all QC metrics are exported in the standardised [mzQC](https://hupo-psi.github.io/mzQC/) format for interoperability.
- **Visualisation** - interactive Plotly plots (barplots, TIC overlays, ion maps, PCA plots, instrument header plots, ...) and a summary table are created over all MS runs.

The identification, the mzQC export and the visualisation can be skipped individually. The mzQC export and visualisation can also be run on the HDF5 files of an earlier run, without repeating the metrics extraction (see [Usage](usage.md#running-only-parts-of-the-pipeline)).

## Documentation

- [Usage](usage.md) - how to run the pipeline, samplesheet format, and all parameters.
- [Parameters](parameters.md) - full parameter reference generated from the pipeline schema.
- [Output](output.md) - description of all output files and directories.
- [Contributing](CONTRIBUTING.md) - how to contribute to the pipeline.
