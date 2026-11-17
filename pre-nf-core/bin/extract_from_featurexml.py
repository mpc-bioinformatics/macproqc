#!/usr/bin/env python

import os
import argparse
import base64
import zipfile
import pyopenms
import h5py
from collections import defaultdict

import hdf5_functions as mzhdf5


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-featurexml", help="FeatureXML with already annotated identifications")
    parser.add_argument("-out_hdf5", help="The Output statistics HDF5")
    parser.add_argument("-report_up_to_charge", type=int, default=5, help="Upper limit of range to be reported in a csv table for the charge")

    return parser.parse_args()


if __name__ == "__main__":
    args = argparse_setup()

    # Count features
    features = pyopenms.FeatureMap()
    pyopenms.FeatureXMLFile().load(args.featurexml, features)

    total_num_features = 0
    num_features_charge = defaultdict(lambda: 0)
    total_num_ident_features = 0
    num_ident_features_charge = defaultdict(lambda: 0) # Identified by charge state
    for f in features:
        # Get data to be able to count
        charge = f.getCharge()
        idents_of_f = f.getPeptideIdentifications()
        is_identified = len(idents_of_f) > 0

        # Count accordingly
        # Count features with charge
        total_num_features += 1
        num_features_charge[charge] += 1

        # Count identified, if any
        if is_identified:
            total_num_ident_features += 1
            num_ident_features_charge[charge] += 1

    zipfile.ZipFile("featurexml.zip", mode="w", compresslevel=9).write(args.featurexml, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9, arcname=args.featurexml.split(os.sep)[-1])
    with open("featurexml.zip", "rb") as fb:
        feature_str_bs64 = base64.b64encode(fb.read())

    with h5py.File(args.out_hdf5, "w") as out_h5:
        mzhdf5.add_entry_to_hdf5(
            f=out_h5,
            qc_acc="MS:4000102",
            qc_short_name="nr_features",
            qc_name="number of detected quantification data points",
            qc_description=(
                "The number of data points detected for quantification purposes within the run. "
                "These data points may be for example XIC profiles, isotopic pattern areas, or reporter ions (see MS:1001805). "
                "The used type should be noted in the metadata or analysis methods section of the recording file for the respective run."
            ),
            value=total_num_features,
            value_shape=(1,),
            value_type="uint64",
            unit_accession="UO:0000189",
            unit_name="count unit",
        )

        mzhdf5.add_entry_to_hdf5(
            f=out_h5,
            qc_acc="MS:4000103",
            qc_short_name="nr_ident_features",
            qc_name="number of identified quantification data points",
            qc_description=(
                "The number of identified data points for quantification purposes within the run after user defined acceptance criteria are applied. "
                "These data points may be for example XIC profiles, isotopic pattern areas, or reporter ions (see MS:1001805). "
                "The used type should be noted in the metadata or analysis methods section of the recording file for the respective run. "
                "In case of multiple acceptance criteria (FDR) available in proteomics, PSM-level FDR should be used for better comparability."
            ),
            value=total_num_ident_features,
            value_shape=(1,),
            value_type="uint64",
            unit_accession="UO:0000189",
            unit_name="count unit",
        )

        report_up_to_charge = args.report_up_to_charge
        charges_features_more = 0
        for key in num_features_charge.keys():
            if key > report_up_to_charge:
                charges_features_more += num_features_charge[key]

        charges_features_frac_list = [
            num_features_charge[i] / total_num_features
            for i in range(1, report_up_to_charge + 1)
        ] + [charges_features_more / total_num_features]

        mzhdf5.add_table_to_hdf5(
            f=out_h5,
            qc_acc="MS:4000207",
            qc_short_name="features_charges",
            qc_name="detected quantification data points charges fractions",
            qc_description=(
                "The fraction of all data points detected for quantification purposes within the run for each specified charge state. "
                "The fractions [0,1] are given in the 'fraction' column, corresponding charges in the 'charge state' column. "
                "The highest charge state is to be interpreted as that charge state or higher."
            ),
            column_names=["MS:1000041 ! charge state", "UO:0000191 ! fraction"],
            column_data=[
                [int(i) for i in range(1, int(report_up_to_charge) + 2)],
                charges_features_frac_list,
            ],
            column_types=["uint64", "float64"],
        )

        charges_features_id_more = 0
        for key in num_ident_features_charge.keys():
            if key > report_up_to_charge:
                charges_features_id_more += num_ident_features_charge[key]

        charges_features_id_frac_list = [
            num_ident_features_charge[i] / total_num_ident_features
            for i in range(1, report_up_to_charge + 1)
        ] + [charges_features_id_more / total_num_ident_features]

        mzhdf5.add_table_to_hdf5(
            f=out_h5,
            qc_acc="MS:4000208",
            qc_short_name="ident_features_charge",
            qc_name="identified quantification data points charges fractions",
            qc_description=(
                "The fraction of all data points detected and identified for quantification purposes within the run for each specified charge state. "
                "The fractions [0,1] are given in the 'fraction' column, corresponding charges in the 'charge state' column. "
                "The highest charge state is to be interpreted as that charge state or higher."
            ),
            column_names=["MS:1000041 ! charge state", "UO:0000191 ! fraction"],
            column_data=[
                [int(i) for i in range(1, int(report_up_to_charge) + 2)],
                charges_features_id_frac_list,
            ],
            column_types=["uint64", "float64"],
        )
