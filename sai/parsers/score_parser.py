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
import re
from sai.parsers.argument_validation import positive_int
from sai.parsers.argument_validation import existed_file
from sai.parsers.argument_validation import between_zero_and_one
from sai.parsers.argument_validation import validate_stat_type
from sai.sai import score
from sai.utils.utils import parse_ind_file


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
        - anc_alleles : str
            Path to the ancestral allele file.
        - w : float
            Allele frequency threshold for the reference group.
        - y : list of float
            List of allele frequency thresholds for each source population. Its length must match `num_src`.
        - output : str
            Path to the output file for storing results.
        - stat_type: str
            Specifies the type of statistic to compute.

    Raises
    ------
    ValueError
        If the length of `args.y` does not match the expected number of source populations (`args.num_src`),
        or if other input parameters do not meet expected conditions.
    """
    src_samples = parse_ind_file(args.src)
    num_src = len(src_samples.keys())
    if len(args.y) != num_src:
        raise ValueError(
            f"The length of y ({len(args.y)}) does not match the number of source populations ({num_src}) found in {args.src}."
        )

    score(
        vcf_file=args.vcf,
        chr_name=args.chr_name,
        ref_ind_file=args.ref,
        tgt_ind_file=args.tgt,
        src_ind_file=args.src,
        win_len=args.win_len,
        win_step=args.win_step,
        num_src=num_src,
        anc_allele_file=args.anc_alleles,
        w=args.w,
        y=args.y,
        output_file=args.output,
        stat_type=args.stat,
        num_workers=1,
    )


def _parse_y_thresholds(value: str) -> tuple[str, float]:
    """
    Parses the --y parameter value to extract an operator and a numerical threshold.

    This function ensures that the input is correctly formatted as one of the following:
    - `=X`  (equality condition)
    - `>X`  (greater than condition)
    - `<X`  (less than condition)
    - `>=X` (greater than or equal to condition)
    - `<=X` (less than or equal to condition)

    The numerical value `X` must be within the range [0, 1].

    Parameters
    ----------
    value : str
        A string representing the allele frequency threshold condition, e.g., "=0.7", ">0.8", "<=0.2".

    Returns
    -------
    tuple[str, float]
        A tuple containing:
        - A string representing the comparison operator (`=`, `<`, `>`, `<=`, `>=`).
        - A float representing the threshold value.

    Raises
    ------
    argparse.ArgumentTypeError
        If the input format is invalid or the numerical threshold is outside the range [0, 1].
    """
    match = re.match(r"^(=|<|>|<=|>=)(\d*\.?\d+)$", value)
    if not match:
        raise argparse.ArgumentTypeError(
            f"Invalid format for --y: {value}. Must be in the form =X, >X, <X, >=X, or <=X "
            f"(e.g., =0.7, >0.8, <0.1, >=0.5, <=0.2)."
        )

    operator, num_str = match.groups()
    num = float(num_str)

    if not (0 <= num <= 1):
        raise argparse.ArgumentTypeError(
            f"Value for --y must be between 0 and 1, got {num}."
        )

    return operator, num


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
        "--ref",
        type=existed_file,
        required=True,
        help="Path to the file with reference population identifiers.",
    )
    parser.add_argument(
        "--tgt",
        type=existed_file,
        required=True,
        help="Path to the file with target population identifiers.",
    )
    parser.add_argument(
        "--src",
        type=existed_file,
        required=True,
        help="Path to the file with source population identifiers.",
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
        "--w",
        type=between_zero_and_one,
        default=0.01,
        help="Frequency threshold for variants in the reference population; only variants with frequencies below this threshold are included in the analysis. Default: 0.01.",
    )
    parser.add_argument(
        "--y",
        type=_parse_y_thresholds,
        nargs="+",
        default=[("=", 1.0)],
        help="List of allele frequency conditions for the source populations. "
        "Each value must be in the form =X, >X, <X, >=X, or <=X "
        "(e.g., =0.7, >0.8, <0.1, >=0.5, <=0.2). "
        "The number of values must match the number of source populations in the file specified by `--src`; "
        "the order of the allele frequency conditions should also correspond to the order of source populations in that file. Default: =1",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output file path for saving results.",
    )
    parser.add_argument(
        "--stat",
        type=validate_stat_type,
        required=True,
        help="Type of statistic to compute: UXX or QXX, where XX is a percentage-like index indicating a threshold in the target population. For example, `U50` means the allele frequency is greater than 0.5, and `Q95` means the allele frequency is greater than or equal to the 95th percentile among sites meeting the specified conditions.",
    )
    parser.set_defaults(runner=_run_score)
