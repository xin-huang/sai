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


import argparse
import sys
from sai.parsers.argument_validation import positive_int
from sai.parsers.argument_validation import positive_number
from sai.parsers.argument_validation import existed_file
from sai.parsers.argument_validation import between_zero_and_one
from sai.sai import outlier


def _run_outlier(args: argparse.Namespace) -> None:
    """
    Runs the outlier detection process based on command-line arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments containing input score file,
        output directory, output prefix, and quantile threshold.
    """
    # Call the outlier function with parsed arguments
    outlier(
        score_file=args.score,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        quantile=args.quantile,
    )


def add_outlier_parser(subparsers: argparse.ArgumentParser) -> None:
    """
    Initializes and configures the command-line interface parser
    for the outlier subcommand.

    Parameters
    ----------
    subparsers : argparse.ArgumentParser
        A command-line interface parser to be configured.
    """
    parser = subparsers.add_parser(
        "outlier", help="Detect and output outlier rows based on quantile thresholds."
    )
    parser.add_argument(
        "--score",
        type=existed_file,
        required=True,
        help="Path to the input score file.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Directory to save the output files.",
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        required=True,
        help="Prefix for the output filenames.",
    )
    parser.add_argument(
        "--quantile",
        type=between_zero_and_one,
        default=0.99,
        help="Quantile threshold for outlier detection, between 0 and 1.",
    )
    parser.set_defaults(runner=_run_outlier)
