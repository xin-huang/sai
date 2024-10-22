# Copyright 2024 Xin Huang
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


import h5py
import multiprocessing
import numpy as np
import pandas as pd


def write_h5(
    file_name: str, group_name: str, data_dict: dict, lock: multiprocessing.Lock
) -> None:
    """
    Writes the data dictionary to an HDF5 file.

    Parameters
    ----------
    file_name : str
        Path to the HDF5 file.
    group_name : str
        Name of the group in the HDF5 file where the data will be stored.
    data_dict : dict
        Dictionary containing the data to be written to the HDF5 file.
    lock : multiprocessing.Lock
        Lock for synchronizing multiprocessing operations.

    """

    converted_dict = {}
    for key, value in data_dict.items():
        if isinstance(value, list):
            if isinstance(value[0], str):  # If the list contains strings
                converted_dict[key] = np.array(
                    value, dtype="S"
                )  # Convert to bytes for HDF5
            else:
                converted_dict[key] = np.array(value)
        else:
            converted_dict[key] = value

    with lock:
        with h5py.File(file_name, "a") as f:
            group = f.create_group(group_name)
            for key, value in converted_dict.items():
                if np.isscalar(value):
                    group.create_dataset(key, data=value)
                else:
                    group.create_dataset(key, data=value, compression="lzf")


def write_tsv(file_name: str, data_dict: dict, lock: multiprocessing.Lock) -> None:
    """
    Writes the data dictionary to a TSV file.

    Parameters
    ----------
    file_name : str
        Path to the TSV file.
    data_dict : dict
        Dictionary containing the data to be written to the TSV file.
    lock : multiprocessing.Lock
        Lock for synchronizing multiprocessing operations.

    """

    converted_dict = {}
    for key, value in data_dict.items():
        if isinstance(value, np.ndarray):
            array_list = value.tolist()
            converted_dict[key] = array_list
        else:
            converted_dict[key] = value

    df = pd.DataFrame([converted_dict])

    with lock:
        with open(file_name, "a") as f:
            df.to_csv(f, header=False, index=False, sep="\t")
