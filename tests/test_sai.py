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
import pandas as pd
import sai.stats
from sai.sai import score, outlier


@pytest.fixture
def example_data(tmp_path):
    # Define example file paths
    pytest.example_vcf = "tests/data/example.vcf"
    pytest.example_config = "tests/data/test_sai.config.yaml"

    # Create a temporary output file path for the score function
    temp_output_file = tmp_path / "output.tsv"

    return {
        "vcf_file": pytest.example_vcf,
        "output_file": str(temp_output_file),
        "output_dir": tmp_path,
        "config": pytest.example_config,
    }


def test_score(example_data):
    # Run score function and capture output
    score(
        vcf_file=example_data["vcf_file"],
        chr_name="21",
        win_len=6666,
        win_step=6666,
        anc_allele_file=None,
        output_file=example_data["output_file"],
        config=example_data["config"],
        num_workers=1,
    )

    # Read the generated output file and validate contents
    df = pd.read_csv(example_data["output_file"], sep="\t")

    col_name = [col for col in df.columns if col.startswith("Q")][0]

    assert df[col_name].iloc[0] == 0.9, "Unexpected value in 'Q' column"


def test_score_mixed_ploidy(example_data):
    score(
        vcf_file="tests/data/test.mixed.ploidy.data.vcf",
        chr_name="21",
        win_len=50000,
        win_step=50000,
        anc_allele_file=None,
        output_file=example_data["output_file"],
        config="tests/data/test_mixed_ploidy.config.yaml",
        num_workers=1,
    )

    df = pd.read_csv(example_data["output_file"], sep="\t")

    assert df["U"].iloc[0] == 0, "Unexpected value in 'U' column"
    assert df["U"].iloc[1] == 1, "Unexpected value in 'U' column"


def test_outlier(example_data):
    output_prefix = f"{example_data['output_dir']}/outliers"
    outlier(
        score_file="tests/data/test.q.scores",
        output_prefix=output_prefix,
        quantile=0.25,
    )

    df = pd.read_csv(f"{output_prefix}.Q.0.25.outliers.tsv", sep="\t")
    assert df["Q"].iloc[0] == 0.7

    outlier(
        score_file="tests/data/test.q.scores",
        output_prefix=output_prefix,
        quantile=0.75,
    )

    df = pd.read_csv(f"{output_prefix}.Q.0.75.outliers.tsv", sep="\t")
    assert df["Q"].iloc[0] == 1.0
    assert df["Q"].iloc[1] == 1.0
