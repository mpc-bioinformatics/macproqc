#!/usr/bin/env python

import sys
import argparse
from alphatims.bruker import TimsTOF
import json


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d_folder", help="Bruker .d-Folder of raw spectra")
    parser.add_argument("-in_json", help="In json query file (similar as to TRFP)")
    parser.add_argument("-out_json", help="Output results file (similar as to TRFP)")
    return parser.parse_args()


if __name__ == "__main__":
    args = argparse_setup()

    # Load queries
    with open(args.in_json, "r") as injson:
        queries = json.load(injson)

    # Open raw spectra
    raw_data = TimsTOF(args.d_folder)

    # Generate header and output format
    result_json = dict()
    result_json["OutputMeta"] = dict(
        base64=False,
        timeunit="minutes",
    )
    result_json["Content"] = []  # Add results as provided by the query

    for entry in queries:
        if entry["tolerance_unit"] != "ppm":
            print("Not supported tolerance unit!")
            sys.exit(1)
        delta = (entry["mz"] / 1000000) * entry["tolerance"]
        min_mz = entry["mz"] - delta
        max_mz = entry["mz"] + delta

        rt_start = entry["rt_start"]*60
        rt_end = entry["rt_end"]*60

        # Indexing as follows: raw_data[RT_in_Secs, Scan_Index, Prec_Index, MZ, Intensity]
        result = dict(
            Meta=dict(MzStart=min_mz, MzEnd=max_mz, RtStart=entry["rt_start"], RtEnd=entry["rt_end"], Comment=entry["comment"]),
            RetentionTimes=list(raw_data[rt_start:rt_end, :, :, min_mz:max_mz, :]["rt_values"]/60),
            Intensities=list(raw_data[rt_start:rt_end, :, :, min_mz:max_mz, :]["intensity_values"])
        )

        result_json["Content"].append(result)

    # Write to file
    with open(args.out_json, "w") as ojs: 
        json.dump(result_json, ojs, indent=2)
