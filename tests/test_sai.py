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


import os
import pytest
import pandas as pd
import tempfile
from sai.sai import score, outlier, plot


@pytest.fixture
def example_data(tmp_path):
    # Define example file paths
    pytest.example_vcf = "./tests/data/example.vcf"
    pytest.example_ref_ind_list = "./tests/data/example.ref.ind.list"
    pytest.example_tgt_ind_list = "./tests/data/example.tgt.ind.list"
    pytest.example_src_ind_list = "./tests/data/example.src.ind.list"

    # Create a temporary output file path for the score function
    temp_output_file = tmp_path / "output.tsv"

    return {
        "vcf_file": pytest.example_vcf,
        "ref_ind_file": pytest.example_ref_ind_list,
        "tgt_ind_file": pytest.example_tgt_ind_list,
        "src_ind_file": pytest.example_src_ind_list,
        "output_file": str(temp_output_file),
        "output_dir": tmp_path,
    }


def test_score(example_data, capfd):
    # Run score function and capture output
    score(
        vcf_file=example_data["vcf_file"],
        chr_name="21",
        ref_ind_file=example_data["ref_ind_file"],
        tgt_ind_file=example_data["tgt_ind_file"],
        src_ind_file=example_data["src_ind_file"],
        win_len=6666,
        win_step=6666,
        num_src=1,
        anc_allele_file=None,
        w=0.3,
        y=[("=", 1)],
        output_file=example_data["output_file"],
        stat_type="Q95",
        num_workers=1,
    )

    # Capture the printed output
    out, err = capfd.readouterr()

    # Assert expected output in stdout
    assert "Found 3 variants with missing genotypes, removing them ..." in out

    # Read the generated output file and validate contents
    df = pd.read_csv(example_data["output_file"], sep="\t")

    col_name = [col for col in df.columns if col.startswith("Q95(")][0]

    assert df[col_name].iloc[0] == 0.9, "Unexpected value in 'Q95' column"


def test_outlier(example_data):
    score(
        vcf_file=example_data["vcf_file"],
        chr_name="21",
        ref_ind_file=example_data["ref_ind_file"],
        tgt_ind_file=example_data["tgt_ind_file"],
        src_ind_file=example_data["src_ind_file"],
        win_len=6666,
        win_step=6666,
        num_src=1,
        anc_allele_file=None,
        w=0.3,
        y=[("=", 1)],
        output_file=example_data["output_file"],
        stat_type="Q95",
        num_workers=1,
    )

    outliers_file = example_data["output_dir"] / "outliers.tsv"

    # Run the outlier function
    with pytest.warns(UserWarning, match="only one unique value"):
        outlier(
            score_file=str(example_data["output_file"]),
            output=str(outliers_file),
            quantile=0.95,
        )

    assert outliers_file.exists()

    q_outliers_file = example_data["output_dir"] / "q_outliers.tsv"

    outlier(
        score_file="tests/data/test.q.scores",
        output=str(q_outliers_file),
        quantile=0.25,
    )

    df = pd.read_csv(str(q_outliers_file), sep="\t")
    assert df["Q95"].iloc[0] == 0.7

    outlier(
        score_file="tests/data/test.q.scores",
        output=str(q_outliers_file),
        quantile=0.75,
    )

    df = pd.read_csv(str(q_outliers_file), sep="\t")
    assert df["Q95"].iloc[0] == 1.0
    assert df["Q95"].iloc[1] == 1.0


def test_plot():
    """Test if the plot function correctly generates an output file."""
    with tempfile.NamedTemporaryFile(delete=False, suffix=".png") as tmp_output:
        output_file = tmp_output.name

    # Call the plot function
    plot(
        u_file="tests/data/test.u.outliers.tsv",
        q_file="tests/data/test.q.outliers.tsv",
        output=output_file,
        xlabel="Q Values",
        ylabel="U Values",
        title="Test Plot",
    )

    # Check if the output file is created
    assert os.path.exists(output_file), "Plot output file was not created."

    # Clean up temporary files
    os.remove(output_file)
