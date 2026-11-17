#!/usr/bin/env python
"""
Helper-functions for the handling of quality metrics in HDF5 files
"""

from typing import List, Any

import numpy as np
import h5py


def add_entry_to_hdf5(
    f: h5py.File,
    qc_acc: str,
    qc_short_name: str,
    qc_name: str,
    qc_description: str,
    value,
    value_shape: tuple,
    value_type: str,
    unit_accession: str = None,
    unit_name: str = None,
) -> None:
    """Adds a single value entry into the hdf5 file"""

    key = f"{qc_acc} ! {qc_short_name}"  # like in OBO for terms, the key is "ACCESSION ! SHORT_NAME"

    if value_type in ("str", h5py.string_dtype()):
        ds = f.create_dataset(
            key, shape=value_shape, dtype=h5py.string_dtype(), compression="gzip"
        )
        ds[:] = value
    else:
        f.create_dataset(key, value_shape, dtype=value_type, compression="gzip")
        if not any([x == 0 for x in value_shape]):
            # Check if any dimension is 0. If so, we do not write data in it (zero lengthed result).
            f[key].write_direct(np.array(value, dtype=value_type))

    f[key].attrs["qc_short_name"] = qc_short_name
    f[key].attrs["qc_name"] = qc_name
    f[key].attrs["qc_description"] = qc_description
    f[key].attrs["unit_accession"] = unit_accession
    f[key].attrs["unit_name"] = unit_name


def add_table_to_hdf5(
    f: h5py.File,
    qc_acc: str,
    qc_short_name: str,
    qc_name: str,
    qc_description: str,
    column_names: List[str],
    column_data: List[List[Any]],
    column_types: List[str],
):
    """Adds a table with an arbitrary number of columns to the HDF5 file."""

    key = f"{qc_acc} ! {qc_short_name}"  # like in OBO for terms, the key is "ACCESSION ! SHORT_NAME"

    table_group = f.create_group(key)
    table_group.attrs["qc_short_name"] = qc_short_name
    table_group.attrs["qc_name"] = qc_name
    table_group.attrs["qc_description"] = qc_description
    table_group.attrs["column_order"] = "|".join(column_names)

    for n, d, t in zip(column_names, column_data, column_types):
        if t in ("str", h5py.string_dtype()):
            ds = table_group.create_dataset(
                n, shape=len(d), dtype=h5py.string_dtype(), compression="gzip"
            )
            ds[:] = d
        else:
            table_group.create_dataset(n, (len(d),), dtype=t, compression="gzip")
            table_group[n].write_direct(np.array(d, dtype=t))
