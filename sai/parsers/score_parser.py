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
        - ref : str
            Path to the reference group individual file.
        - tgt : str
            Path to the target group individual file.
        - src : str
            Path to the source population individual file.
        - win_len : int
            Length of each analysis window.
        - win_step : int
            Step size for moving the window along the sequence.
        - num_src : int
            Number of source populations. The length of `args.y` should match `num_src`.
        - anc_allele_file : str
            Path to the ancestral allele file.
        - ploidy : int
            Ploidy level of the organism.
        - is_phased : bool
            Boolean indicating whether the input data is phased.
        - w : float
            Allele frequency threshold for the reference group.
        - x : float
            Allele frequency threshold for the target group.
        - y : list of float
            List of allele frequency thresholds for each source population. Its length must match `num_src`.
        - output : str
            Path to the output file for storing results.
        - q : float
            Quantile to compute for allele frequencies in the target population.
        - workers : int
            Number of workers to use for parallel processing.

    Raises
    ------
    ValueError
        If the length of `args.y` does not match the expected number of source populations (`args.num_src`),
        or if other input parameters do not meet expected conditions.
    """
    if len(args.y) != args.num_src:
        raise ValueError(
            f"The length of y ({len(args.y)}) does not match the expected number of source populations (num_src = {args.num_src})."
        )

    score(
        vcf_file=args.vcf,
        chr_name=args.chr_name,
        ref_ind_file=args.ref,
        tgt_ind_file=args.tgt,
        src_ind_file=args.src,
        win_len=args.win_len,
        win_step=args.win_step,
        num_src=args.num_src,
        anc_allele_file=args.anc_allele_file,
        ploidy=args.ploidy,
        is_phased=args.is_phased,
        w=args.w,
        x=args.x,
        y=args.y,
        output_file=args.output,
        quantile=args.q,
        num_workers=args.workers,
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
        "score", help="run the score command based on specified parameters."
    )
    parser.add_argument(
        "--vcf",
        type=existed_file,
        required=True,
        help="path to the VCF file containing variant data.",
    )
    parser.add_argument(
        "--chr-name",
        dest="chr_name",
        type=str,
        required=True,
        help="chromosome name to analyze from the VCF file.",
    )
    parser.add_argument(
        "--ref",
        type=existed_file,
        required=True,
        help="path to the file with reference population identifiers.",
    )
    parser.add_argument(
        "--tgt",
        type=existed_file,
        required=True,
        help="path to the file with target population identifiers.",
    )
    parser.add_argument(
        "--src",
        type=existed_file,
        required=True,
        help="path to the file with source population identifiers.",
    )
    parser.add_argument(
        "--win-len",
        dest="win_len",
        type=positive_int,
        default=50000,
        help="length of each genomic window in base pairs. Default is 50,000.",
    )
    parser.add_argument(
        "--win-step",
        dest="win_step",
        type=positive_int,
        default=10000,
        help="step size in base pairs between consecutive windows. Default is 10,000.",
    )
    parser.add_argument(
        "--num-src",
        dest="num_src",
        type=positive_int,
        default=1,
        help="number of source populations to include. Default is 1.",
    )
    parser.add_argument(
        "--anc-allele-file",
        dest="anc_allele_file",
        type=existed_file,
        default=None,
        help="path to the file with ancestral allele information. Default is None.",
    )
    parser.add_argument(
        "--ploidy",
        type=positive_int,
        default=2,
        help="ploidy level of the genome. Default is 2.",
    )
    parser.add_argument(
        "--phased",
        dest="is_phased",
        action="store_true",
        help="specify if genotype data is phased. Default is False.",
    )
    parser.add_argument(
        "--w",
        type=positive_number,
        required=True,
        help="frequency threshold for variants in the reference population; only variants with frequencies below this threshold are included in the analysis.",
    )
    parser.add_argument(
        "--x",
        type=positive_number,
        required=True,
        help="frequency threshold for variants in the target population; only variants with frequencies exceeding this threshold are included in the analysis.",
    )
    parser.add_argument(
        "--y",
        type=between_zero_and_one,
        nargs="+",
        required=True,
        help="list of frequency thresholds for variants in the source populations; only variants with frequencies matching these thresholds are included in the analysis.",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="output file path for saving processed feature vectors.",
    )
    parser.add_argument(
        "--q",
        type=between_zero_and_one,
        default=0.95,
        help="quantile threshold for calculating the q statistic. Default is 0.95.",
    )
    parser.add_argument(
        "--workers",
        type=positive_int,
        default=1,
        help="number of parallel workers for processing. Default is 1.",
    )
    parser.set_defaults(runner=_run_score)
