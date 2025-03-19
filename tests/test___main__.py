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


from unittest.mock import patch, MagicMock
from sai.__main__ import main


@patch("sai.__main__._sai_cli_parser")  # Mock _sai_cli_parser to control its output
@patch(
    "sai.__main__._set_sigpipe_handler"
)  # Mock _set_sigpipe_handler as it doesn’t need testing
def test_main(mock_set_sigpipe_handler, mock_sai_cli_parser):
    # Mock parser and its return values
    mock_parser = MagicMock()
    mock_args = MagicMock()
    mock_args.runner = MagicMock()

    # Configure _sai_cli_parser to return the mock parser
    mock_sai_cli_parser.return_value = mock_parser
    # Configure the mock parser to return mock_args when parse_args is called
    mock_parser.parse_args.return_value = mock_args

    # Call the main function with a test argument list
    test_args = ["score", "--vcf", "tests/data/example.vcf", "--chr-name", "chr1"]
    main(test_args)

    # Check if _set_sigpipe_handler was called
    mock_set_sigpipe_handler.assert_called_once()

    # Verify parse_args was called with test_args
    mock_parser.parse_args.assert_called_once_with(test_args)

    # Ensure runner was called with the parsed arguments
    mock_args.runner.assert_called_once_with(mock_args)
