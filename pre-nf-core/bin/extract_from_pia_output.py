#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import io
import h5py

import hdf5_functions as mzhdf5


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pia_peptides", help="peptides.txt from PIA output")
    parser.add_argument("--pia_proteins", help="Proteins.mzTab from PIA output")
    parser.add_argument("--pia_PSMs", help="PSM.mzTab from PIA output")
    parser.add_argument("--out_hdf5", help="Output HDF5 with statistics")
    parser.add_argument("--store_all_infos", help="Store all PSMs, peptides and proteins in the HDF5 file", default=False, type=bool)

    return parser.parse_args()


def run_pia_extraction():
    args = argparse_setup()

    # Open HDF5 file in write mode
    with h5py.File(args.out_hdf5, 'w') as out_h5:
        parse_protein_infos(args.pia_proteins, out_h5)
        parse_peptide_infos(args.pia_peptides, out_h5)
        parse_psm_infos(args.pia_PSMs, out_h5, args.store_all_infos)


def parse_peptide_infos(pia_peptide_csv: str, out_hdf5: h5py.File) -> int:
    # needs the PIA peptides.csv output file and counts PEPTIDE occurences
    # peptides are already filtered for FDR and decoys by the PIA export
    with open(pia_peptide_csv) as f:
        text = "\n".join([line for line in f if line.startswith("PEPTIDE") or line.startswith("COLS_PEPTIDE") or line.startswith('"COLS_PEPTIDE"')])

    if text:
        peptide_df = pd.read_csv(io.StringIO(text), sep=",")
        number_filtered_peptides = peptide_df.shape[0]
    else:
        number_filtered_peptides = 0
    
    mzhdf5.add_entry_to_hdf5(out_hdf5,
                      "MS:1003250",
                      "nr_peptides",
                      "count of identified peptidoforms",
                      "The number of peptidoforms that pass the threshold to be considered identified with sufficient confidence.", 
                      number_filtered_peptides,
                      (1,),
                      "uint64", 
                      "UO:0000189",
                      "count unit")


def parse_protein_infos(pia_proteins_mztab: str, out_hdf5: h5py.File):
    # read in proteins
    with open(pia_proteins_mztab) as f:
        text = "\n".join([line for line in f if line.startswith("PRH") or line.startswith("PRT")])

    if text:
        protein_df = pd.read_csv(io.StringIO(text), sep="\t")

        # proteins are already FDR and decoy filtered by PIA settings
        number_proteins = protein_df.shape[0]

        ungrouped_proteins = list()
        ungrouped_proteins.extend(protein_df['accession'].values.tolist())
        ungrouped_proteins.extend(protein_df[pd.isna(protein_df['ambiguity_members']) == False]['ambiguity_members'].values.tolist())
        number_ungrouped_proteins = len(set(ungrouped_proteins))
    else:
        number_proteins = 0
        number_ungrouped_proteins = 0

    mzhdf5.add_entry_to_hdf5(
        out_hdf5,
        "MS:1003327",
        "nr_protein_groups",
        "number of identified protein groups",
        "The number of protein groups that pass the threshold to be considered identified with sufficient confidence.",
        number_proteins,
        (1,),
        "uint64",
        "UO:0000189",
        "count unit",
    )

    mzhdf5.add_entry_to_hdf5(
        out_hdf5,
        "MS:4000214",
        "nr_accessions",
        "number of all identified accessions in all ambiguity groups",
        (
            "The number of accessions in identified protein ambiguity groups that have been identified. "
            "This is the number of accessions in the groups, which were counted in 'MS:1002404 ! count of identified proteins', which hence must be greater or equal to this number."
        ),
        number_ungrouped_proteins,
        (1,),
        "uint64",
        "UO:0000189",
        "count unit",
    )


def parse_psm_infos(pia_psm_mztab: str, out_hdf5: h5py.File, store_all_infos: bool = False):
    # get header for PSM FDR Score
    with open(pia_psm_mztab) as f:
        text = "\n".join([line for line in f if line.startswith("MTD")])

    if text: 
        mzTab_header = pd.read_csv(io.StringIO(text), sep="\t", header=None)
        psm_score_header = mzTab_header[mzTab_header[2].str.contains('MS:1002355')].iloc[-1][1] # this is the PSM-level FDR score
        psm_score_header = psm_score_header[4:]     # remove "psm_" at beginning

        ## read in PSMs
        with open(pia_psm_mztab) as f:
            text = "\n".join([line for line in f if line.startswith("PS")])

        psm_df = pd.read_csv(io.StringIO(text), sep="\t")

        # remove decoys and filter for important PSM columns
        psm_df = psm_df.loc[psm_df['opt_global_cv_MS:1002217_decoy_peptide'] == 0]
        psm_df = psm_df.loc[:,["PSM_ID", "sequence", "accession", "unique", "retention_time", "charge", "opt_global_missed_cleavages", "modifications", "exp_mass_to_charge", "calc_mass_to_charge", "spectra_ref", psm_score_header]]

        # group the accessions
        gbseries = psm_df.groupby(by=['PSM_ID'])['accession']
        psm_df['accession'] = psm_df["PSM_ID"].map(gbseries.apply(",".join))
        psm_df = psm_df.drop_duplicates()

        unfiltered_psms = psm_df.copy()
        unfiltered_psms.rename(columns={psm_score_header: 'psm_level_fdrscore'}, inplace=True)

        # filter FDR <= 0.01    # TODO: parameterize?
        psm_df = psm_df.loc[psm_df[psm_score_header] <= 0.01]

        # Calculate ppm error
        exp_calc_diff = psm_df["exp_mass_to_charge"] - psm_df["calc_mass_to_charge"]
        exp_calc_diff_removed_isotopes = exp_calc_diff - (exp_calc_diff.round()).astype(int) # Remove Isotopes, since calc_mass expects none
        ppm_error_df = (exp_calc_diff_removed_isotopes * 1000000) / psm_df["calc_mass_to_charge"]
        ppm_error = ppm_error_df.to_list()

        # now get all counts as QC metrics
        PSM_count = psm_df.shape[0]

        nr_psms = psm_df.shape[0] if psm_df.shape[0] != 0 else 1
        charge_fraction_1    = psm_df[psm_df['charge'] == 1]['PSM_ID'].count() / nr_psms
        charge_fraction_2    = psm_df[psm_df['charge'] == 2]['PSM_ID'].count() / nr_psms
        charge_fraction_3    = psm_df[psm_df['charge'] == 3]['PSM_ID'].count() / nr_psms
        charge_fraction_4    = psm_df[psm_df['charge'] == 4]['PSM_ID'].count() / nr_psms
        charge_fraction_5    = psm_df[psm_df['charge'] == 5]['PSM_ID'].count() / nr_psms
        charge_fraction_more = psm_df[psm_df['charge'] > 5]['PSM_ID'].count() / nr_psms
        charge_fractions = [charge_fraction_1, charge_fraction_2, charge_fraction_3, charge_fraction_4, charge_fraction_5, charge_fraction_more]

        miss_count_0    = psm_df[psm_df['opt_global_missed_cleavages'] == 0]['PSM_ID'].count()
        miss_count_1    = psm_df[psm_df['opt_global_missed_cleavages'] == 1]['PSM_ID'].count()
        miss_count_2    = psm_df[psm_df['opt_global_missed_cleavages'] == 2]['PSM_ID'].count()
        miss_count_3    = psm_df[psm_df['opt_global_missed_cleavages'] == 3]['PSM_ID'].count()
        miss_count_more = psm_df[psm_df['opt_global_missed_cleavages'] > 3]['PSM_ID'].count()
        miss_counts = np.array([miss_count_0, miss_count_1, miss_count_2, miss_count_3, miss_count_more], dtype=np.float64)

    else:
        PSM_count = 0
        charge_fractions = [0, 0, 0, 0, 0, 0]
        miss_counts = miss_counts = np.array([0, 0, 0, 0, 0], dtype=np.float64)
        ppm_error = [np.nan]

    mzhdf5.add_entry_to_hdf5(
        f=out_hdf5,
        qc_acc="MS:1003251",
        qc_short_name="nr_PSMs",
        qc_name="count of identified spectra",
        qc_description="The number of spectra that pass the threshold to be considered identified with sufficient confidence.",
        value=PSM_count,
        value_shape=(1,),
        value_type="uint64",
        unit_accession="UO:0000189",
        unit_name="count unit",
    )

    mzhdf5.add_table_to_hdf5(
        f=out_hdf5,
        qc_acc="MS:4000209",
        qc_short_name="PSM_charge_fractions",
        qc_name="peptide spectrum matches charges fractions",
        qc_description=(
            "The fraction of filtered peptide spectrum matches (PSMs) within the run for each specified charge state. "
            "The fractions [0,1] are given in the 'fraction' column, corresponding charges in the 'charge state' column. "
            "The highest charge state is to be interpreted as that charge state or higher. "
            "The numbers here are recorded after any filtering for PSM level quality (e.g. FDR filtering)."
        ),
        column_names=["MS:1000041 ! charge state", "UO:0000191 ! fraction"],
        column_data=[[1, 2, 3, 4, 5, 6], charge_fractions],
        column_types=["uint16", "float64"],
    )

    miss_fractions = miss_counts / miss_counts.sum() if miss_counts.sum() > 0 else miss_counts
    mzhdf5.add_table_to_hdf5(
        f=out_hdf5,
        qc_acc="MS:4000215",
        qc_short_name="PSM_missed_cleavages_fractions",
        qc_name="missed cleavages fractions",
        qc_description=(
            "The fraction of identified peptides with corresponding number of missed cleavages after user-defined acceptance criteria are applied. "
            "The number of missed cleavages per peptide is given in the 'number of missed cleavages' column, the respective fraction of such peptides identified in the 'fraction' column. "
            "The highest 'missed cleavages' row is to be interpreted as that number of missed cleavages or higher."
        ),
        column_names=[
            "MS:1003044 ! number of missed cleavages",
            "UO:0000191 ! fraction",
        ],
        column_data=[[0, 1, 2, 3, 4], miss_fractions],
        column_types=["uint16", "float64"],
    )

    mzhdf5.add_entry_to_hdf5(
        f=out_hdf5,
        qc_acc="MS:4000206",
        qc_short_name="filtered_psms_ppm_error_quartiles",
        qc_name="precursor ppm deviation distribution",
        qc_description=("The quantiles of the distribution of observed precursor mass accuracies (MS:4000072) [in ppm] of identified MS2 spectra after user-defined acceptance criteria (FDR) are applied. "
            "E.g. one value triplet represents the quartiles Q1, Q2, Q3."),
        value=np.percentile(ppm_error, [25, 50, 75]),   # calculate 1st, 2nd and 3rd quartile of ppm deviations
        value_shape=(3,),
        value_type="float64",
        unit_accession="UO:0000169",
        unit_name="parts per million",
    )

    mzhdf5.add_entry_to_hdf5(
        f=out_hdf5,
        qc_acc="MS:4000178",
        qc_short_name="filtered_psms_ppm_error_mean",
        qc_name="precursor ppm deviation mean",
        qc_description=("The mean of the distribution of observed precursor mass accuracies (MS:4000072) [in ppm] of identified MS2 spectra after user-defined acceptance criteria (FDR) are applied."),
        value=np.mean(ppm_error),
        value_shape=(1,),
        value_type="float64",
        unit_accession="UO:0000169",
        unit_name="parts per million",
    )

    mzhdf5.add_entry_to_hdf5(
        f=out_hdf5,
        qc_acc="MS:4000179",
        qc_short_name="filtered_psms_ppm_error_sigma",
        qc_name="precursor ppm deviation sigma",
        qc_description=("The standard deviation of the distribution of observed precursor mass accuracies (MS:4000072) [in ppm] of identified MS2 spectra after user-defined acceptance criteria (FDR) are applied."),
        value=np.std(ppm_error),
        value_shape=(1,),
        value_type="float64",
        unit_accession="UO:0000169",
        unit_name="parts per million",
    )

    if store_all_infos and (unfiltered_psms is not None):
        # TODO: fix in PIA
        unfiltered_psms["modifications"] = unfiltered_psms["modifications"].astype("str")

        column_names = unfiltered_psms.columns.tolist()
        psms_data = [unfiltered_psms[x].astype("str") for x in column_names]
        column_types = unfiltered_psms.dtypes.astype("str").tolist()

        mzhdf5.add_table_to_hdf5(
            f=out_hdf5,
            qc_acc="LOCAL:identified_psms_table",
            qc_short_name="identified_psms_table",
            qc_name="table of all unfiltered identified PSMs",
            qc_description="All identified PSM as reported by PIA",
            column_names=column_names,
            column_data=psms_data,
            column_types=column_types,
        )

if __name__ == "__main__":
    run_pia_extraction()
