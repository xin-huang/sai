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
from pydantic import ValidationError
from sai.configs import PloidyConfig


def test_valid_config():
    config = PloidyConfig(
        {
            "ref": {"popA": 2},
            "tgt": {"popB": 2},
            "src": {"popC": 2, "popD": 4},
            "outgroup": {"popE": 1},
        }
    )
    assert config.get_ploidy("ref", "popA") == 2
    assert config.get_ploidy("src", "popD") == 4


def test_missing_required_keys():
    with pytest.raises(ValidationError, match="Missing required ploidy keys"):
        PloidyConfig(
            {
                "ref": {"popA": 2},
                "src": {"popC": 2},
            }
        )


def test_extra_keys():
    with pytest.raises(ValidationError, match="Unsupported ploidy keys"):
        PloidyConfig(
            {
                "ref": {"popA": 2},
                "tgt": {"popB": 2},
                "src": {"popC": 2},
                "badkey": {"popD": 3},
            }
        )


def test_invalid_ploidy_values():
    with pytest.raises(ValidationError, match="must be a positive integer"):
        PloidyConfig(
            {
                "ref": {"popA": 2},
                "tgt": {"popB": 0},
                "src": {"popC": 1},
                "outgroup": {"popD": 1},
            }
        )

    with pytest.raises(ValidationError, match="must be a positive integer"):
        PloidyConfig(
            {
                "ref": {"popA": 2},
                "tgt": {"popB": 2},
                "src": {"popC": -2, "popD": -3},
                "outgroup": {"popE": 1},
            }
        )


def test_get_ploidy_key_error():
    config = PloidyConfig(
        {
            "ref": {"popA": 2},
            "tgt": {"popB": 2},
            "src": {"popC": 2},
            "outgroup": {"popD": 1},
        }
    )
    with pytest.raises(KeyError):
        config.get_ploidy("ghost", "popE")
