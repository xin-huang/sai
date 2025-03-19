# Copyright 2025 Xin Huang
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


from multiprocessing import Pool
from typing import Any
from sai.utils.generators import DataGenerator
from sai.utils.preprocessors import DataPreprocessor


def mp_worker(params: tuple[DataPreprocessor, dict]) -> Any:
    """
    Executes the `run` method of the `DataPreprocessor` with provided parameters.

    Parameters
    ----------
    params : tuple of (DataPreprocessor, dict)
        A tuple containing an instance of `DataPreprocessor` and a dictionary of parameters.

    Returns
    -------
    Any
        The result of `data_processor.run(**param_dict)`.
    """
    data_processor, param_dict = params
    return data_processor.run(**param_dict)


def mp_pool(
    data_processor: DataPreprocessor,
    data_generator: DataGenerator,
    nprocess: int,
) -> None:
    """
    Distributes data processing tasks across multiple processes.

    Parameters
    ----------
    data_processor : DataPreprocessor
        An instance of `DataPreprocessor` responsible for processing data.
    data_generator : DataGenerator
        A generator that yields parameter dictionaries for processing.
    nprocess : int
        The number of worker processes to use.

    Returns
    -------
    None
        The processed results are handled by `data_processor.process_items()`.
    """
    tasks: list[tuple[DataPreprocessor, dict]] = [
        (data_processor, params) for params in data_generator.get()
    ]
    with Pool(processes=nprocess) as pool:
        results = pool.map(mp_worker, tasks)

    data_processor.process_items(results)
