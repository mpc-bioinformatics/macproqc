#!/usr/bin/env python
 

import sys
import argparse

import pyopenms
import pandas as pd


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-mztab", help="The mzTab-Identification file")
    parser.add_argument("-out_idxml", help="The Output idxml")

    return parser.parse_args()


if __name__ == "__main__":
    args = argparse_setup()

    # Taken and modified from https://github.com/OpenMS/OpenMS/issues/4417
    run_name = 'unknown'
    peptide_ids = []
    skiplines = 0
    with open(args.mztab) as f_in:
        for line in f_in:
            if line.split('\t', 1)[0] != 'PSH':
                if 'ms_run[1]-location' in line:
                    run_name = line.split('\t')[2]
            else: 
                break
            skiplines += 1

    if skiplines != 0:
        psms = pd.read_csv(args.mztab, sep='\t', header=skiplines-1)
        for _, psm in psms.iterrows():
            peptide_id = pyopenms.PeptideIdentification()
            peptide_id.setRT(psm['retention_time'])
            peptide_id.setMZ(psm['exp_mass_to_charge'])
            peptide_id.setScoreType('q-value')
            peptide_id.setHigherScoreBetter(False) # This was the suggested change from the comment
            # peptide_id.setIdentifier(psm['spectra_ref'])
            peptide_id.setIdentifier(run_name)
            peptide_hit = pyopenms.PeptideHit()
            peptide_hit.setScore(psm['search_engine_score[1]'])
            peptide_hit.setRank(1)
            peptide_hit.setCharge(psm['charge'])
            peptide_hit.setSequence(pyopenms.AASequence.fromString(psm['sequence']))
            peptide_id.setHits([peptide_hit])
            peptide_ids.append(peptide_id)

    protein_id = pyopenms.ProteinIdentification()
    protein_id.setIdentifier(run_name)
    pyopenms.IdXMLFile().store(args.out_idxml, [protein_id], peptide_ids)

