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
from sai.parsers.score_parser import add_score_parser
from sai.parsers.score_parser import _parse_y_thresholds
from sai.parsers.score_parser import _run_score


@pytest.fixture
def parser():
    # Initialize the argument parser with a subparser for the 'score' command
    main_parser = argparse.ArgumentParser()
    subparsers = main_parser.add_subparsers(dest="command")
    add_score_parser(subparsers)
    return main_parser


def test_add_score_parser(parser):
    # Simulate command-line arguments to parse
    args = parser.parse_args(
        [
            "score",
            "--vcf",
            "tests/data/example.vcf",
            "--chr-name",
            "chr1",
            "--ref",
            "tests/data/example.ref.ind.list",
            "--tgt",
            "tests/data/example.tgt.ind.list",
            "--src",
            "tests/data/example.src.ind.list",
            "--win-len",
            "50000",
            "--win-step",
            "10000",
            "--num-src",
            "2",
            "--w",
            "0.3",
            "--x",
            "0.5",
            "--y",
            "=0.1",
            "<0.2",
            "--output",
            "output/results.tsv",
            "--stat",
            "Q95",
        ]
    )

    # Validate parsed arguments
    assert args.command == "score"
    assert args.vcf == "tests/data/example.vcf"
    assert args.chr_name == "chr1"
    assert args.ref == "tests/data/example.ref.ind.list"
    assert args.tgt == "tests/data/example.tgt.ind.list"
    assert args.src == "tests/data/example.src.ind.list"
    assert args.win_len == 50000
    assert args.win_step == 10000
    assert args.num_src == 2
    assert args.anc_alleles is None
    assert args.w == 0.3
    assert args.x == 0.5
    assert args.y == [("=", 0.1), ("<", 0.2)]
    assert args.output == "output/results.tsv"
    assert args.stat == "Q95"


def test_add_score_parser_with_invalid_src():
    # Create an ArgumentParser instance
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    add_score_parser(subparsers)  # Add the `score` subcommand to the parser

    # Simulate command-line input
    args1 = parser.parse_args(
        [
            "score",
            "--vcf",
            "tests/data/example.vcf",
            "--chr-name",
            "chr1",
            "--ref",
            "tests/data/example.ref.ind.list",
            "--tgt",
            "tests/data/example.tgt.ind.list",
            "--src",
            "tests/data/test.src.duplicated.ind.list",  # This file contains only 1 population
            "--win-len",
            "50000",
            "--win-step",
            "10000",
            "--num-src",
            "2",  # Expecting 2 populations
            "--w",
            "0.3",
            "--x",
            "0.5",
            "--y",
            "=0.1",
            "<0.2",
            "--output",
            "output/results.tsv",
            "--stat",
            "U",
        ]
    )

    args2 = parser.parse_args(
        [
            "score",
            "--vcf",
            "tests/data/example.vcf",
            "--chr-name",
            "chr1",
            "--ref",
            "tests/data/example.ref.ind.list",
            "--tgt",
            "tests/data/example.tgt.ind.list",
            "--src",
            "tests/data/example.src.ind.list",
            "--win-len",
            "50000",
            "--win-step",
            "10000",
            "--num-src",
            "2",  # Expecting 2 populations
            "--w",
            "0.3",
            "--x",
            "0.5",
            "--y",
            "=0.1",
            "--output",
            "output/results.tsv",
            "--stat",
            "U",
        ]
    )

    # Trigger _run_score to validate the arguments and raise ValueError
    with pytest.raises(
        ValueError, match=r"The number of populations in the file .* does not match .*"
    ):
        _run_score(args1)

    with pytest.raises(ValueError, match=r"The length of y .* does not match .*"):
        _run_score(args2)


def test_parse_y_thresholds_valid():
    """Test valid cases for _parse_y_thresholds."""
    assert _parse_y_thresholds("=0.7") == ("=", 0.7)
    assert _parse_y_thresholds(">0.8") == (">", 0.8)
    assert _parse_y_thresholds("<0.1") == ("<", 0.1)
    assert _parse_y_thresholds(">=0.5") == (">=", 0.5)
    assert _parse_y_thresholds("<=0.2") == ("<=", 0.2)
    assert _parse_y_thresholds("=1.0") == ("=", 1.0)
    assert _parse_y_thresholds("=0.0") == ("=", 0.0)


def test_parse_y_thresholds_invalid_format():
    """Test cases where input format is invalid."""
    with pytest.raises(argparse.ArgumentTypeError):
        _parse_y_thresholds("0.7")  # Missing operator
    with pytest.raises(argparse.ArgumentTypeError):
        _parse_y_thresholds("==0.7")  # Double equals
    with pytest.raises(argparse.ArgumentTypeError):
        _parse_y_thresholds("!0.5")  # Invalid operator
    with pytest.raises(argparse.ArgumentTypeError):
        _parse_y_thresholds("0.5>")  # Reversed order
    with pytest.raises(argparse.ArgumentTypeError):
        _parse_y_thresholds("")  # Empty string
    with pytest.raises(argparse.ArgumentTypeError):
        _parse_y_thresholds("=<0.3")  # Invalid operator sequence


def test_parse_y_thresholds_out_of_range():
    """Test cases where the numerical value is out of range [0,1]."""
    with pytest.raises(argparse.ArgumentTypeError):
        _parse_y_thresholds("=1.1")  # Out of range
    with pytest.raises(argparse.ArgumentTypeError):
        _parse_y_thresholds(">1.5")
    with pytest.raises(argparse.ArgumentTypeError):
        _parse_y_thresholds("<-0.1")  # Negative value
    with pytest.raises(argparse.ArgumentTypeError):
        _parse_y_thresholds("=-0.2")
