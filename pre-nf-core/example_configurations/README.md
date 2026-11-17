# Files in this folder

## spike-ins.csv

This is a custom format for this QC-Workflow. Which gives us the information which PSMs in the identification have/can be found to also extract the XICs from.

Column-Explanation:

`name`: The actual accession how we could find this entry in an PSM. *NOTE*: We also have to update the FASTA-file to this peptide with the chosen accession.

`sequence`: The sequence which is should be also present in the FASTA. *NOTE*: "old-isa" does not have a sequence and ist not added into the FASTA, hence we can never find identifications for it!

`mz`: The mz, which will be used for the XIC extraction for this entry. We use 10 ppm as a tolerance here

`RT`: The retention time, where we choose the range: (RT / 60) +-3 for XIC-extraction. *Note*: In case of a found PSM, we also extract the XIC of the found PSMs RT. This will be then annotated in a delta-shift, calculated as "found_psms_RT - original_RT"

`mz-tol`: the m/z tolerance in ppm for the spiked in peptide

`rt-tol`: the m/z tolerance in seconds for the spike in
