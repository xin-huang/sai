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


import pytest
import argparse
from sai.parsers.score_parser import add_score_parser


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
            "--anc-allele-file",
            "tests/data/test.anc.allele.bed",
            "--ploidy",
            "2",
            "--phased",
            "--w",
            "0.3",
            "--x",
            "0.5",
            "--y",
            "0.1",
            "0.2",
            "--output",
            "output/results.tsv",
            "--q",
            "0.95",
            "--workers",
            "4",
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
    assert args.anc_allele_file == "tests/data/test.anc.allele.bed"
    assert args.ploidy == 2
    assert args.is_phased is True
    assert args.w == 0.3
    assert args.x == 0.5
    assert args.y == [0.1, 0.2]
    assert args.output == "output/results.tsv"
    assert args.q == 0.95
    assert args.workers == 4
