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
from sai.parsers.score_parser import add_score_parser
from sai.parsers.outlier_parser import add_outlier_parser
from sai.parsers.plot_parser import add_plot_parser


def _set_sigpipe_handler() -> None:
    """
    Sets the signal handler for SIGPIPE signals on POSIX systems.

    """
    import os
    import signal

    if os.name == "posix":
        # Set signal handler for SIGPIPE to quietly kill the program.
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def _sai_cli_parser() -> argparse.ArgumentParser:
    """
    Initializes and configures the command-line interface parser
    for sai.

    Returns
    -------
    top_parser : argparse.ArgumentParser
        A configured command-line interface parser.
    """
    top_parser = argparse.ArgumentParser()
    subparsers = top_parser.add_subparsers(dest="subcommand")
    subparsers.required = True

    add_score_parser(subparsers)
    add_outlier_parser(subparsers)
    add_plot_parser(subparsers)

    return top_parser


def main(arg_list: list = None) -> None:
    """
    Main entry for sai.

    Parameters
    ----------
    arg_list : list, optional
        A list containing arguments for sai. Default: None.
    """
    _set_sigpipe_handler()
    parser = _sai_cli_parser()
    args = parser.parse_args(arg_list)
    args.runner(args)
