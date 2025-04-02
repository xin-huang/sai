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
import pytest
from sai.parsers.argument_validation import positive_int
from sai.parsers.argument_validation import positive_number
from sai.parsers.argument_validation import between_zero_and_one
from sai.parsers.argument_validation import existed_file
from sai.parsers.argument_validation import validate_stat_type


def test_positive_int():
    # Valid positive integer
    assert positive_int("5") == 5

    # Not a positive integer (zero)
    with pytest.raises(argparse.ArgumentTypeError, match="0 is not a positive integer"):
        positive_int("0")

    # Negative integer
    with pytest.raises(
        argparse.ArgumentTypeError, match="-1 is not a positive integer"
    ):
        positive_int("-1")

    # Non-integer string
    with pytest.raises(argparse.ArgumentTypeError, match="abc is not a valid integer"):
        positive_int("abc")


def test_positive_number():
    # Valid positive number
    assert positive_number("3.14") == 3.14

    # Not a positive number (zero)
    with pytest.raises(argparse.ArgumentTypeError, match="0 is not a positive number"):
        positive_number("0")

    # Negative number
    with pytest.raises(
        argparse.ArgumentTypeError, match="-2.5 is not a positive number"
    ):
        positive_number("-2.5")

    # Non-numeric string
    with pytest.raises(argparse.ArgumentTypeError, match="xyz is not a valid number"):
        positive_number("xyz")


def test_between_zero_and_one():
    # Values within range
    assert between_zero_and_one("0.5") == 0.5
    assert between_zero_and_one("0") == 0
    assert between_zero_and_one("1") == 1

    # Values out of range
    with pytest.raises(argparse.ArgumentTypeError, match="1.5 is not between 0 and 1"):
        between_zero_and_one("1.5")

    with pytest.raises(argparse.ArgumentTypeError, match="-0.1 is not between 0 and 1"):
        between_zero_and_one("-0.1")

    # Non-numeric string
    with pytest.raises(
        argparse.ArgumentTypeError, match="not_a_number is not a valid number"
    ):
        between_zero_and_one("not_a_number")


def test_existed_file(tmp_path):
    # Create a temporary file for testing
    temp_file = tmp_path / "temp.txt"
    temp_file.write_text("This is a test file.")

    # Validate an existing file path
    assert existed_file(str(temp_file)) == str(temp_file)

    # Validate a non-existent file path
    with pytest.raises(
        argparse.ArgumentTypeError, match="non_existent_file is not found"
    ):
        existed_file("non_existent_file")


def test_valid_inputs():
    assert validate_stat_type("U50") == "U50"
    assert validate_stat_type("Q05") == "Q05"
    assert validate_stat_type("Q95") == "Q95"
    assert validate_stat_type("Q99") == "Q99"


def test_invalid_inputs():
    with pytest.raises(argparse.ArgumentTypeError):
        validate_stat_type("U")

    with pytest.raises(argparse.ArgumentTypeError):
        validate_stat_type("Q")

    with pytest.raises(argparse.ArgumentTypeError):
        validate_stat_type("Q5")

    with pytest.raises(argparse.ArgumentTypeError):
        validate_stat_type("U100")

    with pytest.raises(argparse.ArgumentTypeError):
        validate_stat_type("Q100")

    with pytest.raises(argparse.ArgumentTypeError):
        validate_stat_type("X50")

    with pytest.raises(argparse.ArgumentTypeError):
        validate_stat_type("Qabc")

    with pytest.raises(argparse.ArgumentTypeError):
        validate_stat_type("")
