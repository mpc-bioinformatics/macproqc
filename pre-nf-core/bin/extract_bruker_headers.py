#!/usr/bin/env python

import argparse
import os
import sqlite3
import h5py
import numpy as np
from alphatims.bruker import TimsTOF

import hdf5_functions as mzhdf5


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d_folder", help="FeatureXML with already annotated identifications")
    parser.add_argument("-out_hdf5", help="The Output statistics HDF5")
    parser.add_argument("-headers_to_parse", "-htp", help="The Headers to parse. Can be applied multiple times", action="append", default=[
        "Vacuum_CurrentFore",
        "Vacuum_Extra4thGauge",
        "Vacuum_CurrentHigh",
        "Vacuum_CurrentFunnel",
        "Digitizer_CurrentTemp",
        "TOF_DeviceTempCurrentValue1",
        "TOF_DeviceTempCurrentValue2"
    ])
    parser.add_argument("-frame_headers_to_parse", "-fhtp", help="The Frame Headers to parse. Can be applied multiple times. If a header is not available it will be skipped.", action="append", default=[
        "Pressure",
    ])
    parser.add_argument("-calibrants_to_retrieve", "-calibrants", help="Calibrants which should be retrieved. In the format MZ:Mobility. E.G.: 922.00978:1.1895", action="append", default=[
        "622.0290:0.9913",
        "922.009798:1.1895",
        "1221.990637:1.3820"
    ])
    parser.add_argument("-calibrants_mz_tolerance", "-cal_mz_tol", help="The MZ Tolerance for the calibrants in Th (m/z). Default: 10", type=float, default=10)
    parser.add_argument("-calibrants_mobility_tolerance", "-cal_mob_tol", help="The Mobility Tolerance for the calibrants in 1/K0 (1/K0). Default: 0.1", type=float, default=0.1)

    return parser.parse_args()


def get_calibrant_info(calibrant_mz, calibrant_mobility, mz_tolerance=10, mobility_tolerance=0.1):
    """
    Gets the calibrants and returns, three arrays: RT, MZ and Mobility values.

    This code was inspired by the AlphaTIMS implementation in the following link:

    https://github.com/MannLabs/alphatims/blob/9fad486536a87c1add095fd4bcf03caca479cc13/alphatims/bruker.py#L2179
    """

    calibrant_lower_mz = calibrant_mz - mz_tolerance
    calibrant_upper_mz = calibrant_mz + mz_tolerance
    calibrant_lower_mobility = calibrant_mobility - mobility_tolerance
    calibrant_upper_mobility = calibrant_mobility + mobility_tolerance

    calibrant_values = br_d[
        :,
        calibrant_lower_mobility: calibrant_upper_mobility,
        slice(0,1),
        calibrant_lower_mz: calibrant_upper_mz,
    ]

    # Get the rows, which have the higest intensity for each retention time for the calibrant within a specific mz and mobility tolerance.
    calibrant_values = calibrant_values.loc[calibrant_values[["rt_values", "intensity_values"]].groupby("rt_values").idxmax()["intensity_values"]]

    calibrant_values_rts = np.array(calibrant_values.index)
    calibrant_values_mzs = np.array(calibrant_values["mz_values"])
    calibrant_values_mobilities = np.array(calibrant_values["mobility_values"])

    return calibrant_values_rts, calibrant_values_mzs, calibrant_values_mobilities


if __name__ == "__main__":

    args = argparse_setup()
    con = sqlite3.connect(args.d_folder + os.sep + "analysis.tdf")
    cur = con.cursor()

    # Get all property definitions
    res = cur.execute(
        "SELECT ID, PermanentName from PropertyDefinitions"
    )
    properties = res.fetchall()
    property_names = [x[1] for x in properties]

    print("Available Properties: ")
    for p in property_names:
        print("\t" + p)

    # Filter to only needed ones
    p_index = []
    p_name = []
    p_col_name = []
    for n in args.headers_to_parse:
        try:
            p_index.append(properties[property_names.index(n)][0])
            p_name.append(n)
        except:
            print("WARNING: Property '{}' not found!".format(n))

    # Open HDF5 file in write mode
    with h5py.File(args.out_hdf5, 'w') as out_h5:

        # Extract data for each frame:
        data_dict = dict()
        for idx, name in zip(p_index, p_name):
            res = cur.execute(
                "SELECT Frame, Value from Properties WHERE  Property = {}".format(idx)
            )
            metadata = res.fetchall()

            data_dict[name] = [x[1]if x[1] is not None else np.nan for x in sorted(metadata, key=lambda x: x[0])]

        # Special CASE: Get MS/MS-Type (and additional headers if available from table Frames)
        frame_columns = [x[1] for x in cur.execute("PRAGMA table_info(Frames);").fetchall()]
        headers_to_retrieve = ["Id", "Time", "MsMsType"]  # Standard headers which always will be extracted.
        
        for h in args.frame_headers_to_parse:
            if h in frame_columns:
                if h not in headers_to_retrieve:
                    headers_to_retrieve.append(h)
            else:
                print("WARNING: Frame Header '{}' not found!".format(h))

        res = cur.execute(
            "SELECT " + ", ".join(headers_to_retrieve) + " from Frames"
        )
        frame_data = res.fetchall()
        sorted_frame_data = sorted(frame_data, key=lambda x: x[0])
        
        # Add to final result table
        column_name = list(data_dict.keys()) + headers_to_retrieve
        column_data = [data_dict[x] for x in data_dict.keys()]
        for c in headers_to_retrieve[:]:
            column_data.append(
                [x[headers_to_retrieve.index(c)] for x in sorted_frame_data]
            )
        column_type = ["float64"]*len(column_name)

        mzhdf5.add_table_to_hdf5(
            f=out_h5,
            qc_acc="BRUKER",
            qc_short_name="Extracted_Headers",
            qc_name="The selected extracted Bruker headers.",
            qc_description=(
                "This table can contain various columns, like 'Temperature' and more. "
                "Depending on the input Bruker-file a column may or may not be present in this table."
            ),
            column_names=column_name,
            column_data=column_data,
            column_types=column_type,
        )

        res.close()
        cur.close()
        con.close()

        # SPECIAL Case about PumpPressure
        # NC_Pump_Pressure
        con = sqlite3.connect(args.d_folder + os.sep + "chromatography-data.sqlite")
        cur = con.cursor()
        res = cur.execute(
            "SELECT Description, Id from TraceSources WHERE Description = 'NC_Pump_Pressure'"
        )
        traces = res.fetchall()

        if len(traces) != 0:
            trace_id = traces[0][1]
            res = cur.execute(
                "SELECT Times, Intensities from TraceChunks WHERE Trace = {}".format(trace_id)
            )
            trace_chunks = res.fetchall() 

            # Extract from binary data
            times = []
            data = []
            for t_chunk, d_chunk in trace_chunks:
                for t_off in range(0, len(t_chunk), 8):
                    times.append(
                        int.from_bytes(t_chunk[t_off:t_off+8], byteorder="little")
                    )
                for t_off in range(0, len(d_chunk), 4):
                    data.append(
                        int.from_bytes(d_chunk[t_off:t_off+4], byteorder="little")
                    )
        else:
            times, data = [np.nan], [np.nan]

        # TODO it is not clear which format the timestamps and the data has. All we know is that the data is provided in little endian
        mzhdf5.add_table_to_hdf5(
            f=out_h5,
            qc_acc="MS:4000210",
            qc_short_name="pump_pressure",
            qc_name="vacuum pump pressure",
            qc_description=(
                "The vacuum pump pressure of a run, defined by the retention times and respectively applied pressures. "
                "The values are similar to the ones saved in MS:1000821, but using a tabular representation."
            ),
            column_names=["MS:1000894 ! retention time", "UO:0000109 ! pressure unit"],
            column_data=[times, data],
            column_types=["float64", "float64"],
        )

        # Add Calibrants Info:
        br_d = TimsTOF(args.d_folder)

        try:
            column_name = ["calibrant_mz", "calibrant_mobility", "observed_calibrant_rt", "observed_calibrant_mz", "observed_calibrant_mobility"]
            column_data = [[],[], [], [], []]
            column_type = ["float64", "float64", "float64", "float64", "float64"]

            for calibrant in args.calibrants_to_retrieve:
                mz, mobility = calibrant.split(":", 1)
                mz = float(mz)
                mobility = float(mobility)

                calibrant_rts, calibrant_mzs, calibrant_mobilities = get_calibrant_info(
                    mz,
                    mobility,
                    mz_tolerance=args.calibrants_mz_tolerance,
                    mobility_tolerance=args.calibrants_mobility_tolerance,
                )

                column_data[0] = np.append(column_data[0], [[mz] * len(calibrant_rts)])
                column_data[1] = np.append(
                    column_data[1], [[mobility] * len(calibrant_rts)]
                )
                column_data[2] = np.append(column_data[2], [calibrant_rts])
                column_data[3] = np.append(column_data[3], [calibrant_mzs])
                column_data[4] = np.append(column_data[4], [calibrant_mobilities])

            mzhdf5.add_table_to_hdf5(
                f=out_h5,
                qc_acc="BRUKE_calibrants",
                qc_short_name="bruker_calibrants",
                qc_name="Bruker calibrant information",
                qc_description=(
                    "Extraction of selected calibrants from the Bruker measurement. "
                    "This table contains the following columns: calibrant_mz, calibrant_mobility --> The observed calibrant mz and mobility AND observed_calibrant_rt, observed_calibrant_mz, observed_calibrant_mobility --> The observed calibrant retention time, mz and mobility."
                ),
                column_names=column_name,
                column_data=column_data,
                column_types=column_type,
            )
        except:
            raise ValueError(
                "The calibrant '{}'could not be retrieved. Is it in the correct format?  ('MZ:Mobility')".format(
                    args.calibrants_to_retrieve
                )
            )
