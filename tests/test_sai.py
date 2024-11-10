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
import pandas as pd
from sai.sai import score, outlier


@pytest.fixture
def example_data(tmp_path):
    # Define example file paths
    pytest.example_vcf = "./tests/data/example.vcf"
    pytest.example_ref_ind_list = "./tests/data/example.ref.ind.list"
    pytest.example_tgt_ind_list = "./tests/data/example.tgt.ind.list"
    pytest.example_src_ind_list = "./tests/data/example.src.ind.list"

    # Create a temporary output file path for the score function
    temp_output_file = tmp_path / "output.tsv"

    # Run score function once to create the output file
    score(
        vcf_file=pytest.example_vcf,
        chr_name=21,
        ref_ind_file=pytest.example_ref_ind_list,
        tgt_ind_file=pytest.example_tgt_ind_list,
        src_ind_file=pytest.example_src_ind_list,
        win_len=3333,
        win_step=3333,
        num_src=1,
        anc_allele_file=None,
        ploidy=2,
        is_phased=True,
        w=0.3,
        x=0.5,
        y=[1],
        output_file=str(temp_output_file),
        quantile=0.95,
        num_workers=1,
    )

    return {
        "output_file": temp_output_file,
        "output_dir": tmp_path,
    }


def test_score(example_data):
    # Read the generated output file and validate contents
    df = pd.read_csv(example_data["output_file"], sep="\t")

    assert df['U'].iloc[0] == 3
    assert df['Q95'].iloc[0] == 0.9


def test_outlier(example_data):
    # Define the output directory and prefix for outlier results
    output_dir = example_data["output_dir"] / "outliers"
    output_prefix = "test"

    # Run the outlier function using the score file generated in test_score
    outlier(
        score_file=str(example_data["output_file"]),
        output_dir=str(output_dir),
        output_prefix=output_prefix,
        quantile=0.95,
    )

    # Check if the outlier files are created
    u_outliers_file = output_dir / f"{output_prefix}_U_outliers.tsv"
    q_outliers_file = output_dir / f"{output_prefix}_Q95_outliers.tsv"

    assert u_outliers_file.exists()
    assert q_outliers_file.exists()

    # Optionally, read and check contents of the outlier files
    u_outliers_df = pd.read_csv(u_outliers_file, sep="\t")
    q_outliers_df = pd.read_csv(q_outliers_file, sep="\t")

    assert u_outliers_df.empty, "U outliers file is unexpectedly not empty"
    assert q_outliers_df.empty, "Q outliers file is unexpectedly not empty"
