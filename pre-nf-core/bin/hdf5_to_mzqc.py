#!/usr/bin/env python
import argparse
from datetime import datetime
import os
import re
import json
from typing import Any, Dict, List

import h5py
import pandas as pd

from mzqc import MZQCFile as qc

psi_accession_regex = r"(MS:[A-Z0-9]+).*"  # regular expression for a PSI-MS accession
accession_regex = r"([A-Z0-9]+:[A-Z0-9]+) ! .+"  # regular expression for an accession ONTOLOGY:ACCESSION


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-hdf5", help="Path to the hdf5 file created by McQuaC", required=True)
    parser.add_argument("-mzqc_out", help="Output path for the generated mzQC file", required=True
    )
    return parser.parse_args()


def hdf5_entry_to_mzqc_metric(hdf5_file: h5py.File, key: str) -> qc.QualityMetric:
    """Convert an entry from the hdf5 file to an mzQC element"""

    # split the key (i.e. name of table in hdf) to get the basic accession and short name
    split_key = key.split("!")

    print(f"Splitted key: {split_key} (original key: {key})")

    metric_accession = split_key[0].strip()
    # possible attributes: 'qc_description', 'qc_name', 'qc_short_name', 'unit_accession', 'unit_name'
    metric_name = hdf5_file[key].attrs.get("qc_name", split_key[1].strip())

    metric_unit = None
    unit_accession = hdf5_file[key].attrs.get("unit_accession")
    if (unit_accession is not None) and len(unit_accession) > 0:
        unit_name = hdf5_file[key].attrs.get("unit_name")
        metric_unit = {"unit_accession": unit_accession, "unit_name": unit_name}

    if type(hdf5_file[key]) == h5py.Dataset:
        if hdf5_file[key].shape[0] == 1:
            metric_value = hdf5_file[key][0]
        else:
            metric_value = list(hdf5_file[key])
    elif type(hdf5_file[key]) == h5py.Group:
        # check, if all column names are an accession, defined by the accession_regex
        key_list = list(hdf5_file[key].keys())

        if all(re.fullmatch(accession_regex, key) for key in key_list):
            metric_value = table_from_hdf5_group(hdf5_file[key])
        else:
            metric_value = None

    qm = qc.QualityMetric(
        accession=metric_accession,
        name=metric_name,
        value=metric_value,
        unit=metric_unit,
    )
    return qm


def table_from_hdf5_group(group: h5py.Group) -> Dict[str, List[Any]]:
    """
    Creates an mzQC table from an HDF5 group, assuming that the group contains a table structure.
    Returned is a dictionary with the column names as keys and the values as lists.
    """
    column_names = list(group.keys())

    ret_dict = {}

    for k in column_names:
        dict_key = re.match(accession_regex, k).group(1)
        # Convert bytestrings to normal strings if needed
        dict_val = [
            val.decode("utf-8") if isinstance(val, bytes) else val for val in group[k]
        ]
        ret_dict[dict_key] = dict_val

    return ret_dict


def process_hdf5_to_run_quality(hdf5_path: str, analysis_software: List) -> qc.RunQuality:
    """Process a single HDF5 file and return a RunQuality object."""

    input_filename = os.path.splitext(os.path.basename(hdf5_path))[0]
    fileFormat = qc.CvParameter(accession="MS:1000563", name="Thermo RAW format")

    input_file_start_time_epoch = 0
    quality_metrics = []
    with h5py.File(hdf5_path, "r") as hdf5_file:
        # go through the entries in the hdf5 file and get all the metrics
        for key in hdf5_file.keys():
            if re.fullmatch(psi_accession_regex, key):
                # only try to parse PSI-MS accessions
                qm = hdf5_entry_to_mzqc_metric(hdf5_file, key)

                # only use "valid" metrics for now
                if qm.accession.startswith("MS"):
                    quality_metrics.append(qm)
            elif str(key).startswith("LOCAL:startTime"):
                # get start time of raw file creation (in timestamp)
                input_file_start_time_epoch = hdf5_file[key][0]

    file_properties = []
    if input_file_start_time_epoch > 0:
        readable_time = datetime.fromtimestamp(input_file_start_time_epoch).strftime(
            "%Y-%m-%d %H:%M:%S"
        )
        file_properties.append(
            qc.CvParameter(
                accession="MS:1002435",
                name="data processing start time",
                value=readable_time,
            )
        )
    # here we could add more information about the files, if we had it

    input_file_raw = qc.InputFile(
        name=input_filename,
        location=None,
        fileFormat=fileFormat,
        fileProperties=file_properties,
    )

    meta = qc.MetaDataParameters(
        inputFiles=[input_file_raw], analysisSoftware=analysis_software
    )
    return qc.RunQuality(metadata=meta, qualityMetrics=quality_metrics)


if __name__ == "__main__":
    args = argparse_setup()

    comet_version = "v2024.01.0"

    # we could add information about our search engine
    comet_mzqc = qc.AnalysisSoftware(
        accession="MS:1002251",
        name="Comet",
        description="Comet open-source sequence search engine developed at the University of Washington.",
        version=comet_version,
        uri="https://uwpr.github.io/Comet/",
    )

    cv_ms = qc.ControlledVocabulary(
        name="Proteomics Standards Initiative Mass Spectrometry Ontology",
        version="4.1.197",
        uri="https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo",
    )

    # TODO: and should add information about the feature detection as well
    # anso_nb = qc.AnalysisSoftware(version="0.1.2.3", uri="file:///mylocal/jupyter/host")

    # TODO: and about MCQuaC itself
    # anso_nb = qc.AnalysisSoftware(version="0.1.2.3", uri="file:///mylocal/jupyter/host")

    # create one mzQC file per HDF5 file
    rq = process_hdf5_to_run_quality(args.hdf5, [comet_mzqc])

    mzqc = qc.MzQcFile(
        version="1.0.0",
        creationDate=datetime.now().isoformat(),
        runQualities=[rq],
        setQualities=[],
        controlledVocabularies=[cv_ms],
    )

    with open(args.mzqc_out, "w") as mzqc_file:
        mzqc_file.write(
            json.dumps(json.loads(qc.JsonSerialisable.to_json(mzqc)), indent=2)
        )
