#!/usr/bin/env python

import argparse
from collections import defaultdict
import h5py
import numpy as np

from fisher_py import RawFile
from fisher_py.data import Device, ToleranceUnits
from fisher_py.data.business import (ChromatogramSignal,
                                     ChromatogramTraceSettings,
                                     GenericDataTypes, Scan,
                                     SpectrumPacketType, TraceType)
from fisher_py.data.filter_enums import MsOrderType
from fisher_py.mass_precision_estimator import PrecisionEstimate
from fisher_py.raw_file_reader import RawFileAccess, RawFileReaderAdapter

import hdf5_functions as mzhdf5


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-raw", help="FeatureXML with already annotated identifications")
    parser.add_argument("-out_hdf5", help="The Output statistics HDF5")
    parser.add_argument("-extra_headers_to_parse", "-ehtp", help="The Headers to parse. Can be applied multiple times", action="append", 
        default=[
            "Ion Injection Time (ms)",
            "Number of Lock Masses",
            "Lock Mass #1 (m/z)",
            "Lock Mass #2 (m/z)",
            "Lock Mass #3 (m/z)",
            "LM Search Window (ppm)",
            "LM Search Window (mmu)",
            "Last Locking (sec)",
            "LM m/z-Correction (ppm),LM Correction",
    ])
    parser.add_argument("-tune_headers_to_parse", "-thtp", help="The Headers to parse. Can be applied multiple times. NOTE: These entries are one dimensional.", action="append", 
        default=[
            "Ion Transfer Tube Temperature (+ or +-)",
            "Ion Transfer Tube Temperature (-)",
            "Vaporizer Temp. (+ or +-)",
            "Vaporizer Temp. (-)",
    ])
    parser.add_argument("-log_headers_to_parse", "-lhtp", help="The Headers to parse. Can be applied multiple times", action="append", 
        default=[
            "Ion Transfer Tube Temperature (°C)",
            "Vaporizer Temperature (°C)",
            "Ambient temp. (°C)",
            "Orbitrap block temp. (°C)",
            "Detector temp. (°C)",
            "Ion tr. tube temp. (°C)",
            "Vaporizer temp. (°C)",
            "CPU core temp. (°C)",
            "PCB temp. top (°C)",
            "PCB temp. center (°C)",
            "PCB temp. bottom (°C)",
            "Inner Electrode Temp. (°C)",
            "Outer Electrode 1 Temp. (°C)",
            "Outer Electrode 2 Temp. (°C)",
            "MCB Ambient Temp. (°C)",
            "MCB PCB Temp. (°C)",
            "TX PCB temperature (°C)",
    ])
    return parser.parse_args()


def get_headers_to_parse(headers, headers_from_raw):
    """ Helper function to map which header from which index should be retrieved"""
    statistics_to_retrieve = []
    for idx, h in enumerate(headers_from_raw):
        for hp in headers:
            for hp_part in hp.split(","):
                if h.label.startswith(hp_part):
                    # Retrievable in this RAW-file
                    statistics_to_retrieve.append((idx, hp))
                    break

    return statistics_to_retrieve


if __name__ == "__main__":
    args = argparse_setup()

    data_dict = defaultdict(lambda: list())
    raw_file = RawFileReaderAdapter.file_factory(args.raw)
    raw_file.select_instrument(Device.MS, 1)  # Selecting the MS

    log_statistics_to_retrieve = get_headers_to_parse(
        args.log_headers_to_parse,
        raw_file.get_status_log_header_information()
    )
    tune_statistics_to_retrieve = get_headers_to_parse(
        args.tune_headers_to_parse,
        raw_file.get_tune_data_header_information()
    )
    extra_statistics_to_retrieve = get_headers_to_parse(
        args.extra_headers_to_parse,
        raw_file.get_trailer_extra_header_information()
    )

    # Open HDF5 file in write mode
    with h5py.File(args.out_hdf5, 'w') as out_h5:

        # Retrieve all the information
        first_scan_number = raw_file.run_header_ex.first_spectrum
        last_scan_number = raw_file.run_header_ex.last_spectrum
        for scan in range(first_scan_number, last_scan_number + 1):
            scan_statistics = raw_file.get_scan_stats_for_scan_number(scan)
            start_time_of_scan = scan_statistics.start_time
            data_dict["Scan_StartTime"].append(start_time_of_scan)
            scan_filter = raw_file.get_filter_for_scan_number(scan)
            data_dict["Scan_msLevel"].append(scan_filter.ms_order.value)

            # Get Info of filtered statistics we want to track (log and extra data)
            log_scan_values = raw_file.get_status_log_for_retention_time(start_time_of_scan).values
            for idx, hp in log_statistics_to_retrieve:
                data_dict["LOG_" + hp].append(log_scan_values[idx])
            extra_scan_values = raw_file.get_trailer_extra_information(scan).values
            for idx, hp in extra_statistics_to_retrieve:
                data_dict["EXTRA_" + hp].append(extra_scan_values[idx])

        if len(list(data_dict.keys())) != 0:
            column_name = list(data_dict.keys())
            column_data = [data_dict[x] for x in column_name]
            column_type = ["float64"]*len(column_name)
            mzhdf5.add_table_to_hdf5(
                f=out_h5,
                qc_acc="THERMO",
                qc_short_name="Extracted_Headers",
                qc_name="The selected extracted Thermo headers.",
                qc_description=(
                    "This table can contain various columns, ranging from 'Temperature', 'Lock Masses' and more. "
                    "Depending on the input RAW-file a column may be present in this table."
                ),
                column_names=[x.replace("/", "") for x in column_name],
                column_data=column_data,
                column_types=column_type,
            )

        # Tune data is one dimensional, therefore single values
        tune_scan_values = raw_file.get_tune_data(0).values
        if len(tune_statistics_to_retrieve) != 0:
            tune_dict = dict()
            for idx, hp in tune_statistics_to_retrieve:
                tune_dict["TUNE_" + hp] = tune_scan_values[idx]
            column_name = list(tune_dict.keys())
            column_data = [tune_dict[x] for x in column_name]
            column_type = ["float64"]*len(column_name)

            mzhdf5.add_table_to_hdf5(
                f=out_h5,
                qc_acc="THERMO_LOG",
                qc_short_name="Extracted_Log_Headers",
                qc_name="The selected extracted Thermo log headers.",
                qc_description=(
                    "This table contains various columns of length one, like 'Vaporizer Temperature'. "
                    "Depending on the input RAW-file a column may be present in this table."
                ),
                column_names=[x.replace("/", "") for x in column_name],
                column_data=column_data,
                column_types=column_type,
            )

        # Get all the information from the instruments (in FreeStyle under Devices)
        num_devices = raw_file.get_instrument_count_of_type(Device.Analog)  # Get the number of devices
        settings = ChromatogramTraceSettings(TraceType.Analog1)
        pump_pressure_dict = defaultdict(lambda: list())

        for i in range(1, num_devices+1):
            # Iterate over each device
            raw_file.select_instrument(Device.Analog, i)
            label = raw_file.get_instrument_data().axis_label_y

            # Check if it is the pump_preasure, which we want to extract
            if label in ("Pump_Pressure bar", "NC_Pump_Pressure bar"):
                # This works for QexHF, QeXI and FLI (NC_Pump_Pressure) and EX, EXI and EXII (Pump_Preasure)
                data = raw_file.get_chromatogram_data([settings], 1, -1)
                pump_pressure_dict["pump_pressure_bar_y_axis"] = list(data.intensities_array[0])
                pump_pressure_dict["pump_pressure_bar_x_axis"] = list(data.positions_array[0])

        if "pump_pressure_bar_y_axis" not in pump_pressure_dict:
            # This is currently the case for PROETD and OEI. These are not captured, hence None
            pump_pressure_dict["pump_pressure_bar_x_axis"] = [np.nan]
            pump_pressure_dict["pump_pressure_bar_y_axis"] = [np.nan]

        column_name = ["pump_pressure_x_axis", "pump_pressure_y_axis"]
        column_data = [pump_pressure_dict["pump_pressure_bar_x_axis"], pump_pressure_dict["pump_pressure_bar_y_axis"]]
        column_type = ["float64", "float64"]

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
            column_data=[
                pump_pressure_dict["pump_pressure_bar_x_axis"],
                pump_pressure_dict["pump_pressure_bar_y_axis"],
            ],
            column_types=["float64", "float64"],
        )

    # Close raw-file
    raw_file.dispose()  
