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


import pytest
import numpy as np
from sai.stats.features import calc_u, calc_q
from sai.utils.generators import WindowDataGenerator
from sai.utils.preprocessors import FeatureVectorsPreprocessor


@pytest.fixture
def feature_preprocessor():
    # Create an instance of FeatureVectorsPreprocessor with thresholds and temporary output file
    return FeatureVectorsPreprocessor(
        w=0.3, x=0.5, y=[0.2, 0.4], output_file="test_output.tsv"
    )


def test_run(feature_preprocessor):
    # Define mock input data to test the run method
    chr_name = "21"
    ref_pop = "ref1"
    tgt_pop = "tgt1"
    src_pop_list = ["src1", "src2"]
    start, end = 1000, 2000

    # Create mock genotype data
    ref_gts = np.array([[0, 0, 1], [1, 1, 0], [0, 1, 1]])
    tgt_gts = np.array([[0, 1, 1], [1, 1, 1], [0, 0, 1]])
    src_gts_list = [
        np.array([[0, 0, 0], [1, 0, 0], [1, 1, 1]]),
        np.array([[1, 1, 1], [0, 1, 1], [0, 0, 1]]),
    ]

    # Run the method
    result = feature_preprocessor.run(
        chr_name=chr_name,
        ref_pop=ref_pop,
        tgt_pop=tgt_pop,
        src_pop_list=src_pop_list,
        start=start,
        end=end,
        ref_gts=ref_gts,
        tgt_gts=tgt_gts,
        src_gts_list=src_gts_list,
        ploidy=1,
    )

    # Check that the result contains the expected keys
    assert result["chr_name"] == chr_name
    assert result["start"] == start
    assert result["end"] == end
    assert result["ref_pop"] == ref_pop
    assert result["tgt_pop"] == tgt_pop
    assert result["src_pop_list"] == src_pop_list
    assert "u_statistic" in result
    assert "q_statistic" in result


def test_process_items(feature_preprocessor, tmp_path):
    # Generate a temporary output file path using tmp_path
    temp_output = tmp_path / "test_output.tsv"
    feature_preprocessor.output_file = str(temp_output)

    # Define a mock items dictionary to test process_items method
    items = {
        "chr_name": "21",
        "start": 1000,
        "end": 2000,
        "ref_pop": "ref1",
        "tgt_pop": "tgt1",
        "src_pop_list": ["src1", "src2"],
        "u_statistic": 5,
        "q_statistic": 0.75,
    }

    # Run the process_items method
    feature_preprocessor.process_items(items)

    # Check the output file content
    with open(temp_output, "r") as f:
        lines = f.readlines()
        assert len(lines) == 1  # Ensure only one line is written
        expected_output = "21\t1000\t2000\tref1\ttgt1\tsrc1,src2\t5\t0.75\n"
        assert lines[0] == expected_output


@pytest.fixture
def example_data():
    pytest.example_vcf = "./tests/data/example.vcf"
    pytest.example_ref_ind_list = "./tests/data/example.ref.ind.list"
    pytest.example_tgt_ind_list = "./tests/data/example.tgt.ind.list"
    pytest.example_src_ind_list = "./tests/data/example.src.ind.list"


def test_run_from_file(example_data, tmp_path):
    # Set up the WindowDataGenerator
    generator = WindowDataGenerator(
        vcf_file=pytest.example_vcf,
        chr_name=21,
        ref_ind_file=pytest.example_ref_ind_list,
        tgt_ind_file=pytest.example_tgt_ind_list,
        src_ind_file=pytest.example_src_ind_list,
        win_len=3333,
        win_step=3333,
    )

    # Create a temporary output file path using tmp_path
    temp_output_file = tmp_path / "output.tsv"

    # Initialize the FeatureVectorsPreprocessor with the temporary output file path
    preprocessor = FeatureVectorsPreprocessor(
        w=0.3, x=0.5, y=[1], output_file=str(temp_output_file)
    )

    # Run the generator and preprocessor
    for window_data in generator.get():
        items = preprocessor.run(
            chr_name=window_data["chr_name"],
            ref_pop=window_data["ref_pop"],
            tgt_pop=window_data["tgt_pop"],
            src_pop_list=window_data["src_pop_list"],
            start=window_data["start"],
            end=window_data["end"],
            ref_gts=window_data["ref_gts"],
            tgt_gts=window_data["tgt_gts"],
            src_gts_list=window_data["src_gts_list"],
            ploidy=window_data["ploidy"],
        )
        preprocessor.process_items(items)
        assert items["u_statistic"] == 3
        assert items["q_statistic"] == 0.9

    # Check that the temporary file was created and is not empty
    assert temp_output_file.exists()
    assert temp_output_file.stat().st_size > 0

    # Optionally, open and inspect contents of temp_output_file here if necessary
    with open(temp_output_file, "r") as f:
        lines = f.readlines()
        assert len(lines) > 0  # Ensure some content was written
