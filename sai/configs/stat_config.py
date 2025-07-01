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


SUPPORTED_STATISTICS = [
    "Danc",
    "DD",
    "df",
    "Dplus",
    "fd",
    "U",
    "Q",
]


from pydantic import RootModel, field_validator, ValidationError
from typing import Dict, Literal, List, Optional, Union


class StatConfig(
    RootModel[
        Dict[
            str,
            Dict[str, Dict[str, Union[float, str]]],
        ]
    ]
):
    """
    A class to represent the configuration for various statistics used in the analysis.

    This class manages the configuration of statistical parameters for different
    statistical tests (e.g., "U", "Q"). It validates the range of parameters
    such as `ref`, `tgt`, and `src`, where `ref` and `tgt` are numerical values
    representing frequencies between 0 and 1, and `src` can be a list of strings with
    comparison operators (e.g., "=0.5", ">0.2").
    """

    @field_validator("root")
    def check_valid_stat_types(
        cls, v: Dict[str, Dict[str, Dict[str, Union[float, str]]]]
    ) -> Dict[
        str,
        Dict[str, Dict[str, Union[float, tuple[str, float]]]],
    ]:
        """
        Validates statistics parameters, specifically for U and Q types.

        Parameters
        ----------
        v : Dict[str, Dict[str, Dict[str, Union[float, str]]]]
            A dictionary mapping statistic names (e.g., "U", "Q") to parameter groups ("ref", "tgt", "src"),
            where each group is a mapping of population names to values.

            - Outer dict key: statistic name (e.g., "U", "Q", "fd")
            - Middle dict key: parameter group ("ref", "tgt", or "src")
            - Inner dict key: population name (e.g., "AFR", "CHB")
            - Inner dict value:
                - For "ref" and "tgt": float (frequency between 0 and 1)
                - For "src": string comparator expression (e.g., ">=0.2", "=1")

        Returns
        -------
        Dict[str, Dict[str, Dict[str, Union[float, tuple[str, float]]]]]
            A validated and normalized statistics dictionary.

            - Outer dict key: statistic name (e.g., "U", "Q")
            - Middle dict key: parameter group ("ref", "tgt", "src")
            - Inner dict key: population name (e.g., "AFR", "CHB")
            - Inner dict value:
                - For "ref" and "tgt": float (validated to be between 0 and 1)
                - For "src": tuple (comparator operator, float), e.g., (">=", 0.2)

        Raises
        ------
        ValueError
            If any name of statistics is not supported.
        """
        for stat_name, params in v.items():
            if stat_name not in SUPPORTED_STATISTICS:
                raise ValueError(f"The {stat_name} statistic is not supported.")
            if stat_name in ["U", "Q"]:
                # Validate U and Q statistics parameters
                cls.check_range_for_u_q(stat_name, params)
        return v

    @staticmethod
    def check_range_for_u_q(
        stat_name: str, params: Dict[str, Dict[str, Union[float, str]]]
    ) -> None:
        """
        Validates the parameters for U and Q statistics.
        ref and tgt must be between 0 and 1, and src must contain a valid comparator
        with a frequency value.

        Parameters
        ----------
        stat_name : str
            The name of the statistic (e.g., "U" or "Q").
        params : Dict[str, Dict[str, Union[float, str]]]
            A dictionary containing the parameters for the statistic, such as ref,
            tgt, and src.

        Raises
        ------
        ValueError
            If any of the parameters are outside the valid range or in an incorrect
            format.
        """
        if stat_name in ["U", "Q"]:
            required_keys = {"ref", "tgt", "src"}
            param_keys = set(params.keys())
            if param_keys != required_keys:
                raise ValueError(
                    f"{stat_name} must have exactly the keys: {required_keys}, but got {param_keys}."
                )

        for param, pop_values in params.items():
            if param in ["ref", "tgt"]:
                for pop, value in pop_values.items():
                    num = float(value)
                    if not (0 <= num <= 1):
                        raise ValueError(
                            f"{param}[{pop}] value must be between 0 and 1 for {stat_name}, got {val}."
                        )
            elif param == "src":
                new_src: Dict[str, tuple[str, float]] = {}
                for pop, expr in pop_values.items():
                    if not isinstance(expr, str):
                        raise ValueError(
                            f"{param}[{pop}] value must be a comparator string for {stat_name}."
                        )
                    new_src[pop] = StatConfig.check_comparator(
                        expr, stat_name, f"src[{pop}]"
                    )
                params["src"] = new_src

    @staticmethod
    def check_comparator(value: str, stat_name: str, param: str) -> tuple[str, float]:
        """
        Validates that the src parameter contains a valid comparator (e.g., "=0.5", ">=0.2"),
        and ensure the number is between 0 and 1.

        Parameters
        ----------
        value : str
            The value of the src parameter, which should contain a comparator (e.g., "=0.5").
        stat_name : str
            The name of the statistic (e.g., "U" or "Q").
        param : str
            The parameter name ("src").

        Returns
        -------
        tuple[str, float]
            A tuple containing:
            - A string representing the comparison operator (`=`, `<`, `>`, `<=`, `>=`).
            - A float representing the threshold value.

        Raises
        ------
        ValueError
            If the value does not contain a valid comparator or the number is not in
            the range 0-1.
        """
        valid_comparators = ["<=", ">=", "=", "<", ">"]
        if not any(comp in value for comp in valid_comparators):
            raise ValueError(
                f"{param} for {stat_name} must contain a valid comparator (e.g., '=0.5', '>=0.2')."
            )

        # Extract the numeric value after the comparator
        comparator = next(comp for comp in valid_comparators if comp in value)
        try:
            num = float(value[len(comparator) :])
        except ValueError:
            raise ValueError(
                f"{param} value for {stat_name} must be a valid number after the comparator."
            )

        if not (0 <= num <= 1):
            raise ValueError(
                f"{param} value must be between 0 and 1 for {stat_name}, but got {num}."
            )

        return comparator, num

    def get_parameters(
        self, stat_name: str
    ) -> Optional[Dict[str, Dict[str, Union[float, tuple[str, float]]]]]:
        """
        Retrieves the parameters for a specific statistic.

        Parameters
        ----------
        stat_name : str
            The name of the statistic whose parameters are to be retrieved.

        Returns
        -------
        Optional[Dict[str, Dict[str, Union[float, tuple[str, float]]]]]
            A dictionary containing the parameters for the specified statistic,
            or None if not found.
        """
        return self.root.get(stat_name, None)
