#!/usr/bin/env python

import csv
import json
import argparse
from statistics import mean
from pyteomics import proforma

def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-icsv", help="Input-CSV-File of the SpikeIns")
    parser.add_argument("-iidents", help="Input-Identification-Files in mzTAB format (here the SpikeIns are specially encoded (see FASTA))")
    parser.add_argument("-ojson", help="The JSON-Output-File to be used with the ThermoRawFileParser XIC")
    parser.add_argument("-oidentifications", help="Output mapping the sequences to number of identifications, if a sequence is not listed, there was no (valid) ID")

    return parser.parse_args()

if __name__ == "__main__":

    args = argparse_setup()

    seq_to_data = dict()

    with open(args.icsv, "r") as in_csv_file:
        # Read spike-ins CSV and get header indicies
        csv_in = csv.reader(in_csv_file)
        header = next(csv_in)
        name_idx = header.index("name")
        seq_idx = header.index("sequence")
        mz_idx = header.index("mz")
        rt_idx = header.index("RT")
        mz_tol_idx = header.index("mz-tol")
        rt_tol_idx = header.index("rt-tol")

        # Generate the two output structures
        for l in csv_in:
            # go line-wise through the spike-ins CSV
            
            mz_tol_split = str(l[mz_tol_idx]).split(sep=" ")
            rt_tol = float(l[rt_tol_idx]) / 60
            
            # transform the proforma sequence to a plain sequence (only amino acids, no modifications) to be able to check for identifications, which are only annotated with plain sequences
            proforma_sequence = proforma.parse(str(l[seq_idx]))
            plain_sequence = "".join(aa for aa, mods in proforma_sequence[0])
            
            # map from the sequence to all other data (basically our json)
            seq_to_data[plain_sequence] = {
                "comment": str(l[name_idx]),
                "mz": float(l[mz_idx]),
                "tolerance": float(mz_tol_split[0]),
                "tolerance_unit": str(mz_tol_split[1]),
                "rt_start": (float(l[rt_idx]) / 60) - rt_tol,
                "rt_end": (float(l[rt_idx]) / 60) + rt_tol
            }

    found_rts = dict()
    with open(args.iidents, "r") as in_ident_file:
        rt_idx = 0
        accession_idx = 0
        score_label = ""
        score_idx = 0
        for l in in_ident_file:
            # iterate through mzTab

            # Get Header Line
            if l.startswith("MTD"):
                cols = l.split("\t")
                if cols[1].startswith("psm_search_engine_score") and ("MS:1002257" in cols[2]):
                    # this is the "comet expression" score, remove the "psm_" prefix
                    score_label = cols[1][4:]
            elif l.startswith("PSH"):
                cols = l.split("\t")
                seq_idx = cols.index("sequence") 
                rt_idx = cols.index("retention_time")
                accession_idx = cols.index("accession")
                if len(score_label) > 0:
                    score_idx = cols.index(score_label)

            # Check if we found an psms which actually fits to the SpikeIns and put it into the found RTs
            if l.startswith("PSM"):
                cols = l.split("\t")
                if cols[seq_idx] in seq_to_data.keys():
                    if score_idx < 1 or float(cols[score_idx]) < 0.01:
                        # either cannot filter for score or it is below 0.01
                        if cols[seq_idx] not in found_rts.keys():
                            found_rts[cols[seq_idx]] = []
                        
                        found_rts[cols[seq_idx]].append( float(cols[rt_idx]) )

        # adjust the RTs of found identifications
        for seq, vals in found_rts.items():
            new_rt = mean(set(vals)) / 60

            tol = (seq_to_data[seq]["rt_end"] - seq_to_data[seq]["rt_start"]) / 2

            seq_to_data[seq]["rt_start"] = new_rt - tol
            seq_to_data[seq]["rt_end"] = new_rt + tol
    
    for seq in seq_to_data.keys():
        seq_to_data[seq]["rt_start"] = max(0, seq_to_data[seq]["rt_start"])
    
    # Write Output JSON-File
    with open(args.ojson, "w") as o_json:
        o_json.write(json.dumps([val for val in seq_to_data.values()], indent=4))

    # Write Output identifications file
    with open(args.oidentifications, "w") as output_identifications:
        for seq, vals in found_rts.items():
            nr_ids = len(set(vals))
            output_identifications.write(seq + "," + str(nr_ids) + "\n")
