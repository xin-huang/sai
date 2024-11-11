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


import pytest
import argparse
from sai.parsers.plot_parser import add_plot_parser


@pytest.fixture
def parser():
    # Initialize the argument parser with a subparser for the 'outlier' command
    main_parser = argparse.ArgumentParser()
    subparsers = main_parser.add_subparsers(dest="command")
    add_plot_parser(subparsers)
    return main_parser


def test_add_plot_parser(parser):
    # Simulate command-line arguments to parse
    args = parser.parse_args(
        [
            "plot",
            "--outlier",
            "tests/data/test.outliers.tsv",
            "--output",
            "output.png",
            "--xlabel",
            "Test X Label",
            "--ylabel",
            "Test Y Label",
            "--title",
            "Test Title",
            "--figsize-x",
            "8",
            "--figsize-y",
            "5",
            "--dpi",
            "200",
            "--alpha",
            "0.5",
        ]
    )

    # Check if the attributes are set correctly
    assert args.outlier == "tests/data/test.outliers.tsv"
    assert args.output == "output.png"
    assert args.xlabel == "Test X Label"
    assert args.ylabel == "Test Y Label"
    assert args.title == "Test Title"
    assert args.figsize_x == 8
    assert args.figsize_y == 5
    assert args.dpi == 200
    assert args.alpha == 0.5
