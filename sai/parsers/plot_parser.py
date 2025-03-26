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
from sai.parsers.argument_validation import positive_number
from sai.parsers.argument_validation import existed_file
from sai.sai import plot


def _run_plot(args: argparse.Namespace) -> None:
    """
    Runs the plotting process based on command-line arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments containing input files, output file,
        xlabel, ylabel, title, figsize_x, figsize_y, dpi, alpha,
        marker_size, marker_color, and marker_style.
    """
    plot(
        u_file=args.u_file,
        q_file=args.q_file,
        output=args.output,
        xlabel=args.xlabel,
        ylabel=args.ylabel,
        title=args.title,
        figsize_x=args.figsize_x,
        figsize_y=args.figsize_y,
        dpi=args.dpi,
        alpha=args.alpha,
        marker_size=args.marker_size,
        marker_color=args.marker_color,
        marker_style=args.marker_style,
    )


def add_plot_parser(subparsers: argparse.ArgumentParser) -> None:
    """
    Initializes and configures the command-line interface parser
    for the plot subcommand.

    Parameters
    ----------
    subparsers : argparse.ArgumentParser
        A command-line interface parser to be configured.
    """
    parser = subparsers.add_parser(
        "plot", help="Generate a scatter plot of U vs Q statistics."
    )
    parser.add_argument(
        "--u-file",
        dest="u_file",
        type=existed_file,
        required=True,
        help="Path to the U score/outlier file.",
    )
    parser.add_argument(
        "--q-file",
        dest="q_file",
        type=existed_file,
        required=True,
        help="Path to the Q score/outlier file.",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to save the output plot file. The format depends on the file extension (e.g., `.png`, `.pdf`).",
    )
    parser.add_argument(
        "--xlabel",
        type=str,
        default="Q Statistic",
        help="Label for the X-axis. Default: Q Statistic.",
    )
    parser.add_argument(
        "--ylabel",
        type=str,
        default="U Statistic",
        help="Label for the Y-axis. Default: U Statistic.",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="Scatter Plot of U vs Q",
        help="Title of the plot. Default: Scatter Plot of U vs Q.",
    )
    parser.add_argument(
        "--figsize-x",
        type=positive_number,
        default=6,
        help="Width of the figure (in inches). Default: 6.",
    )
    parser.add_argument(
        "--figsize-y",
        type=positive_number,
        default=6,
        help="Height of the figure (in inches). Default: 6.",
    )
    parser.add_argument(
        "--dpi",
        type=positive_int,
        default=300,
        help="Resolution of the saved plot. Default: 300.",
    )
    parser.add_argument(
        "--alpha",
        type=positive_number,
        default=0.6,
        help="Transparency level of scatter points. Default: 0.6.",
    )
    parser.add_argument(
        "--marker-size",
        dest="marker_size",
        type=positive_number,
        default=20,
        help="Size of the scatter plot markers. See matplotlib.pyplot.scatter. Default: 20.",
    )
    parser.add_argument(
        "--marker-color",
        dest="marker_color",
        type=str,
        default="blue",
        help="Color of the markers. See matplotlib.pyplot.scatter. Default: blue.",
    )
    parser.add_argument(
        "--marker-style",
        dest="marker_style",
        type=str,
        default="o",
        help="Shape of the markers. See matplotlib.pyplot.scatter. Default: o.",
    )
    parser.set_defaults(runner=_run_plot)
