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
from sai.parsers.argument_validation import positive_int
from sai.parsers.argument_validation import existed_file
from sai.sai import score


def _run_score(args: argparse.Namespace) -> None:
    """
    Executes the scoring function with arguments provided via the command line interface.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments that contain the necessary parameters for the scoring function,
        including:

        - vcf : str
            Path to the VCF file containing variant data.
        - chr_name : str
            Name of the chromosome to be analyzed.
        - win_len : int
            Length of each analysis window.
        - win_step : int
            Step size for moving the window along the sequence.
        - anc_alleles : str
            Path to the ancestral allele file.
        - output : str
            Path to the output file for storing results.
        - stat_config: str
            Path to the YAML configuration file specifying the statistics, ploidy levels, and populations to compute.

    Raises
    ------
    ValueError
        If fewer than three ploidy values are provided,
        or if the number of ploidy values for source populations does not match `num_src`.
        or if other input parameters do not meet expected conditions.
    """
    score(
        vcf_file=args.vcf,
        chr_name=args.chr_name,
        win_len=args.win_len,
        win_step=args.win_step,
        anc_allele_file=args.anc_alleles,
        output_file=args.output,
        config=args.config,
        num_workers=1,
    )


def add_score_parser(subparsers: argparse.ArgumentParser) -> None:
    """
    Initializes and configures the command-line interface parser
    for the score subcommand.

    Parameters
    ----------
    subparsers : argparse.ArgumentParser
        A command-line interface parser to be configured.
    """
    parser = subparsers.add_parser(
        "score", help="Run the score command based on specified parameters."
    )
    parser.add_argument(
        "--vcf",
        type=existed_file,
        required=True,
        help="Path to the VCF file containing variant data.",
    )
    parser.add_argument(
        "--chr-name",
        dest="chr_name",
        type=str,
        required=True,
        help="Chromosome name to analyze from the VCF file.",
    )
    parser.add_argument(
        "--win-len",
        dest="win_len",
        type=positive_int,
        default=50000,
        help="Length of each genomic window in base pairs. Default: 50,000.",
    )
    parser.add_argument(
        "--win-step",
        dest="win_step",
        type=positive_int,
        default=10000,
        help="Step size in base pairs between consecutive windows. Default: 10,000.",
    )
    parser.add_argument(
        "--anc-alleles",
        dest="anc_alleles",
        type=existed_file,
        default=None,
        help="Path to the BED file with ancestral allele information. If ancestral allele information is not provided, filtering will be performed for each variant based on whether the allele frequency of any allele (assuming biallelic) meets the specified condition during the calculation of the statistics. Default: None.",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output file path for saving results.",
    )
    parser.add_argument(
        "--config",
        type=existed_file,
        required=True,
        help="Path to the YAML configuration file specifying the statistics to compute, ploidy settings, and population group file paths.",
    )
    parser.set_defaults(runner=_run_score)
