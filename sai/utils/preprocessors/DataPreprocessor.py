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


from typing import Any
from abc import ABC, abstractmethod


class DataPreprocessor(ABC):
    """
    Abstract base class for preprocessing genomic data.

    This class defines a common interface for various data preprocessing operations,
    such as filtering, normalization, and transformation of genomic data. Subclasses
    should implement specific methods to handle data processing tasks, ensuring a
    consistent way to run operations and manage the output of processed data.

    Methods:
    --------
    run(**kwargs) -> Any:
        Execute the core data processing task. Subclasses must define this method
        to carry out specific preprocessing tasks such as filtering, normalization,
        or transformation. This method should return the processed data, which
        will then be handled by the main process to manage further steps or output.

    process_items(items, **kwargs) -> None:
        Handle the output or further processing of data once the `run` method
        has completed. This allows subclasses to define how processed data
        should be managed, such as saving results to a file, database, or converting
        the data to a specific format for future analysis.
    """

    @abstractmethod
    def run(self, **kwargs) -> Any:
        """
        Abstract method to run the preprocessing operations.

        Subclasses must implement this method to perform specific preprocessing
        tasks based on the initialized parameters and any additional keyword
        arguments.

        Parameters:
        -----------
        **kwargs : dict
            Additional keyword arguments that may be required for specific
            preprocessing operations.

        Returns:
        --------
        processed_data : Any
            The result of the preprocessing task, which can be further handled
            by the `process_items` method.
        """
        pass

    @abstractmethod
    def process_items(self, items: Any, **kwargs) -> None:
        """
        Abstract method to handle the output or post-processing of data.

        Subclasses must implement this method to define how the processed data
        should be managed. This could include saving the data to a file,
        transforming it into a new format, or preparing it for the next step
        of analysis.

        Parameters:
        -----------
        items : Any
            The processed data returned by the `run` method, which will be managed
            or output according to the logic defined in this method.

        **kwargs : dict
            Additional keyword arguments that can be used for customizing the
            output process. For example, this may include options like `output_file`
            to specify where the data should be saved or other settings to control
            the output format.
        """
        pass
