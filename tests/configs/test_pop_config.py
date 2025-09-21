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


import tempfile
import os
import pytest
from sai.configs import PopConfig


@pytest.fixture
def temp_files():
    files = {}
    for key in ["ref", "tgt", "src", "outgroup"]:
        f = tempfile.NamedTemporaryFile(delete=False)
        f.write(b"sample1\nsample2\n")
        f.close()
        files[key] = f.name
    yield files
    for path in files.values():
        os.unlink(path)


def test_valid_pop_config_with_outgroup(temp_files):
    config = PopConfig(
        {
            "ref": temp_files["ref"],
            "tgt": temp_files["tgt"],
            "src": temp_files["src"],
            "outgroup": temp_files["outgroup"],
        }
    )
    assert config.get_population("ref") == temp_files["ref"]
    assert config.get_population("outgroup") == temp_files["outgroup"]


def test_valid_pop_config_without_outgroup(temp_files):
    config = PopConfig(
        {
            "ref": temp_files["ref"],
            "tgt": temp_files["tgt"],
            "src": temp_files["src"],
        }
    )
    assert config.get_population("tgt") == temp_files["tgt"]
    assert config.get_population("outgroup") is None


def test_missing_required_key(temp_files):
    with pytest.raises(ValueError, match="Missing required population keys"):
        PopConfig({"ref": temp_files["ref"], "src": temp_files["src"]})


def test_invalid_extra_key(temp_files):
    with pytest.raises(ValueError, match="Unsupported population keys"):
        PopConfig(
            {
                "ref": temp_files["ref"],
                "tgt": temp_files["tgt"],
                "src": temp_files["src"],
                "ghost": temp_files["outgroup"],
            }
        )


def test_file_not_exist():
    with pytest.raises(ValueError, match="does not exist"):
        PopConfig(
            {
                "ref": "/non/existent/path",
                "tgt": "/non/existent/path",
                "src": "/non/existent/path",
            }
        )
