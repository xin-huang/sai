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
from sai.configs.global_config import GlobalConfig
from sai.configs.stat_config import StatConfig
from sai.configs.ploidy_config import PloidyConfig
from sai.configs.pop_config import PopConfig


def test_global_config_valid(tmp_path):
    ref_file = tmp_path / "ref.txt"
    tgt_file = tmp_path / "tgt.txt"
    src_file = tmp_path / "src.txt"

    ref_file.write_text("popA sample1\npopA sample2")
    tgt_file.write_text("popB sample3\npopB sample4")
    src_file.write_text("popC sample5\npopC sample6")

    stat_cfg = StatConfig.model_validate(
        {
            "Q": {
                "ref": {"popA": 0.3},
                "tgt": {"popB": 0.95},
                "src": {"popC": "=1"},
            }
        }
    )

    ploidy_cfg = PloidyConfig.model_validate(
        {
            "ref": {"popA": 2},
            "tgt": {"popB": 2},
            "src": {"popC": 2},
        }
    )

    pop_cfg = PopConfig.model_validate(
        {
            "ref": str(ref_file),
            "tgt": str(tgt_file),
            "src": str(src_file),
        }
    )

    global_cfg = GlobalConfig(
        statistics=stat_cfg,
        ploidies=ploidy_cfg,
        populations=pop_cfg,
    )

    assert global_cfg.statistics.get_parameters("Q")["ref"]["popA"] == 0.3
    assert global_cfg.ploidies.get_ploidy("ref", "popA") == 2
    assert global_cfg.populations.get_population("ref") == str(ref_file)


def test_global_config_invalid_ploidy(tmp_path):
    ref_file = tmp_path / "ref.txt"
    tgt_file = tmp_path / "tgt.txt"
    src_file = tmp_path / "src.txt"

    ref_file.write_text("popA sample1\npopA sample2")
    tgt_file.write_text("popB sample3\npopB sample4")
    src_file.write_text("popC sample5\npopC sample6")

    stat_cfg = StatConfig.model_validate(
        {
            "Q": {
                "ref": {"popA": 0.3},
                "tgt": {"popB": 0.95},
                "src": {"popC": "=1"},
            }
        }
    )

    ploidy_cfg = PloidyConfig.model_validate(
        {
            "ref": {"popA": 2},
            "tgt": {"popD": 2},
            "src": {"popC": 2},
        }
    )

    pop_cfg = PopConfig.model_validate(
        {
            "ref": str(ref_file),
            "tgt": str(tgt_file),
            "src": str(src_file),
        }
    )

    with pytest.raises(
        ValueError,
        match=r"Population 'popB' used in statistics\[Q\]\[tgt\] is not defined in ploidies\[tgt\]",
    ):
        GlobalConfig(
            statistics=stat_cfg,
            ploidies=ploidy_cfg,
            populations=pop_cfg,
        )


def test_global_config_invalid_population(tmp_path):
    ref_file = tmp_path / "ref.txt"
    tgt_file = tmp_path / "tgt.txt"
    src_file = tmp_path / "src.txt"

    ref_file.write_text("popA sample1\npopA sample2")
    tgt_file.write_text("popD sample3\npopD sample4")
    src_file.write_text("popC sample5\npopC sample6")

    stat_cfg = StatConfig.model_validate(
        {
            "Q": {
                "ref": {"popA": 0.3},
                "tgt": {"popB": 0.95},
                "src": {"popC": "=1"},
            }
        }
    )

    ploidy_cfg = PloidyConfig.model_validate(
        {
            "ref": {"popA": 2},
            "tgt": {"popB": 2},
            "src": {"popC": 2},
        }
    )

    pop_cfg = PopConfig.model_validate(
        {
            "ref": str(ref_file),
            "tgt": str(tgt_file),
            "src": str(src_file),
        }
    )

    with pytest.raises(
        ValueError,
        match=r"Population 'popB' used in statistics\[Q\]\[tgt\] is not found in the population file for group 'tgt'",
    ):
        GlobalConfig(
            statistics=stat_cfg,
            ploidies=ploidy_cfg,
            populations=pop_cfg,
        )
