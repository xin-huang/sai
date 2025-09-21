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
import yaml
from sai.configs import StatConfig


def test_stat_config_valid():
    config = StatConfig(
        {
            "U": {
                "ref": {"popA": 0.01},
                "tgt": {"popB": 0.5},
                "src": {"Nea": "=1", "Den": ">0.8"},
            },
            "Q": {
                "ref": {"popA": 0.01},
                "tgt": {"popB": 0.95},
                "src": {"Nea": ">0.2", "Den": "<0.8"},
            },
            "fd": {},  # No parameters
            "df": {},  # No parameters
        }
    )

    # Test valid configuration
    assert config.get_parameters("U") == {
        "ref": {"popA": 0.01},
        "tgt": {"popB": 0.5},
        "src": {"Nea": ("=", 1.0), "Den": (">", 0.8)},
    }

    assert config.get_parameters("Q") == {
        "ref": {"popA": 0.01},
        "tgt": {"popB": 0.95},
        "src": {"Nea": (">", 0.2), "Den": ("<", 0.8)},
    }
    assert config.get_parameters("fd") == {}
    assert config.get_parameters("df") == {}


def test_stat_config_invalid():
    # Test unsupported statistics
    with pytest.raises(ValueError):
        StatConfig({"qq": {}})

    # Test invalid src format (not a list)
    with pytest.raises(ValueError):
        StatConfig(
            {"U": {"ref": "0.01", "tgt": "0.5", "src": "=1"}}  # src should be a list
        )

    # Test out-of-range numeric values for 'ref'
    with pytest.raises(ValueError):
        StatConfig(
            {
                "U": {"ref": -0.1, "tgt": 0.5, "src": ["=1", ">0.8"]},
            }
        )

    # Test out-of-range numeric values for 'tgt'
    with pytest.raises(ValueError):
        StatConfig(
            {
                "U": {"ref": 0.01, "tgt": -100, "src": ["=1", ">0.8"]},
            }
        )

    # Test out-of-range numeric values for 'src'
    with pytest.raises(ValueError):
        StatConfig(
            {
                "U": {"ref": 0.01, "tgt": 0.5, "src": ["=-1", ">0.8"]},
            }
        )

    # Test non-numeric string values for 'ref'
    with pytest.raises(ValueError):
        StatConfig(
            {
                "U": {"ref": "foo", "tgt": 0.5, "src": ["=1", ">0.8"]},
            }
        )

    # Test non-numeric string values for 'tgt'
    with pytest.raises(ValueError):
        StatConfig(
            {
                "U": {"ref": 0.01, "tgt": "foo", "src": ["=1", ">0.8"]},
            }
        )

    # Test non-numeric string values after comparator
    with pytest.raises(ValueError):
        StatConfig(
            {
                "U": {"ref": 0.01, "tgt": 0.5, "src": ["=invalid", ">0.8"]},
            }
        )

    # Test missing ref, tgt, src in U and Q statistics
    with pytest.raises(ValueError):
        StatConfig({"U": {"ref": "0.01", "tgt": "0.5"}})  # src is missing

    with pytest.raises(ValueError):
        StatConfig({"Q": {}})

    # Test invalid src comparator (invalid value)
    with pytest.raises(ValueError):
        StatConfig(
            {
                "Q": {
                    "ref": "0.01",
                    "tgt": "0.95",
                    "src": ["invalid_value"],
                }
            }
        )


def test_stat_config_from_file():
    with open("tests/data/test_config.yaml", "r") as f:
        data = yaml.safe_load(f)

    stat_config = StatConfig(data["statistics"])

    stat_names = list(stat_config.root.keys())
    assert "U" in stat_names
    assert "Q" in stat_names

    u_params = stat_config.get_parameters("U")
    assert u_params == {
        "ref": {"popA": 0.01},
        "tgt": {"popB": 0.5},
        "src": {"Nea": ("=", 1), "Den": (">=", 0.8)},
    }

    q_params = stat_config.get_parameters("Q")
    assert q_params == {
        "ref": {"popA": 0.01},
        "tgt": {"popB": 0.95},
        "src": {"Nea": (">", 0.2), "Den": ("<=", 0.8)},
    }

    assert "ref" in u_params
    assert "tgt" in u_params
    assert "src" in u_params

    with open("tests/data/test_invalid_stat_config.yaml", "r") as f:
        invalid_data = yaml.safe_load(f)

    with pytest.raises(ValueError):
        StatConfig(statistics=invalid_data["statistics"])
