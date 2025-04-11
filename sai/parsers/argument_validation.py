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


import argparse
import os
import re


def positive_int(value: str) -> int:
    """
    Validates if the provided string represents a positive integer.

    Parameters
    ----------
    value : str
        The value to validate.

    Returns
    -------
    int
        The validated positive integer.

    Raises
    ------
    argparse.ArgumentTypeError
        If the value is not a valid integer or positive integer.
    """
    if value is not None:
        try:
            value = int(value)
        except ValueError:
            raise argparse.ArgumentTypeError(f"{value} is not a valid integer")
        if value <= 0:
            raise argparse.ArgumentTypeError(f"{value} is not a positive integer")
    return value


def positive_number(value: str) -> float:
    """
    Validates if the provided string represents a positive number.

    Parameters
    ----------
    value : str
        The value to validate.

    Returns
    -------
    float
        The validated positive number.

    Raises
    ------
    argparse.ArgumentTypeError
        If the value is not a valid number or positive number.
    """
    if value is not None:
        try:
            value = float(value)
        except ValueError:
            raise argparse.ArgumentTypeError(f"{value} is not a valid number")
        if value <= 0:
            raise argparse.ArgumentTypeError(f"{value} is not a positive number")
    return value


def between_zero_and_one(value: str) -> float:
    """
    Validates if the provided string represents a number between 0 and 1 (inclusive).

    Parameters
    ----------
    value : str
        The value to validate.

    Returns
    -------
    float
        The validated number between 0 and 1.

    Raises
    ------
    argparse.ArgumentTypeError
        If the value is not a valid number or is not between 0 and 1.
    """
    if value is not None:
        try:
            value = float(value)
        except ValueError:
            raise argparse.ArgumentTypeError(f"{value} is not a valid number")
        if not (0 <= value <= 1):
            raise argparse.ArgumentTypeError(
                f"{value} is not between 0 and 1 (inclusive)"
            )
    return value


def existed_file(value: str) -> str:
    """
    Validates if the provided string is a path to an existing file.

    Parameters
    ----------
    value : str
        The path to validate.

    Returns
    -------
    str
        The validated file path.

    Raises
    ------
    argparse.ArgumentTypeError
        If the file does not exist.
    """
    if value is not None:
        if not os.path.isfile(value):
            raise argparse.ArgumentTypeError(f"{value} is not found")
    return value


def validate_stat_type(value: str) -> str:
    """
    Validate the input `stat_type`.

    Parameters
    ----------
    value : str
        The statistic type to validate. Must be either:
        - "U" : Compute the U statistic.
        - "QXX" : Compute the Q statistic, where "XX" is a one or two-digit integer
          representing the quantile percentage (e.g., "Q95" for 95th quantile).

    Returns
    -------
    str
        The validated `stat_type`, either "U" or "QXX".

    Raises
    ------
    argparse.ArgumentTypeError
        If the input does not match the expected format ("U" or "QXX").
    """
    if re.fullmatch(
        r"[UQ]\d{2}", value
    ):  # Matches U or Q followed by exactly two digits
        return value
    else:
        raise argparse.ArgumentTypeError(
            f"Invalid --stat-type: {value}. Must be 'UXX' or 'QXX' (e.g., 'U05' for x > 0.05, 'Q95' for quantile = 0.95)."
        )
