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
import sai.stats
from sai.preprocessors import ChunkPreprocessor
from sai.configs import PloidyConfig, StatConfig


@pytest.fixture
def example_data(tmp_path):
    # Define example file paths
    pytest.example_vcf = "tests/data/example.vcf"
    pytest.example_ref_ind_list = "tests/data/example.ref.ind.list"
    pytest.example_tgt_ind_list = "tests/data/example.tgt.ind.list"
    pytest.example_src_ind_list = "tests/data/example.src.ind.list"
    pytest.example_config = "tests/data/example.config.yaml"

    # Create a temporary output file path for the score function
    temp_output_file = tmp_path / "output.tsv"

    return {
        "vcf_file": pytest.example_vcf,
        "ref_ind_file": pytest.example_ref_ind_list,
        "tgt_ind_file": pytest.example_tgt_ind_list,
        "src_ind_file": pytest.example_src_ind_list,
        "config": pytest.example_config,
        "output_file": str(temp_output_file),
        "output_dir": tmp_path,
    }


def test_chunk_preprocessor(example_data):
    with open(example_data["config"], "r") as f:
        config = yaml.safe_load(f)

    stat_config = StatConfig(config["statistics"])
    ploidy_config = PloidyConfig(config["ploidies"])

    preprocessor = ChunkPreprocessor(
        vcf_file=example_data["vcf_file"],
        ref_ind_file=example_data["ref_ind_file"],
        tgt_ind_file=example_data["tgt_ind_file"],
        src_ind_file=example_data["src_ind_file"],
        out_ind_file=None,
        win_len=6666,
        win_step=6666,
        num_src=1,
        anc_allele_file=None,
        output_file=example_data["output_file"],
        stat_config=stat_config,
        ploidy_config=ploidy_config,
    )

    results = preprocessor.run(
        chr_name="21",
        start=0,
        end=6666,
    )

    assert results[0]["Q"] == 0.9
