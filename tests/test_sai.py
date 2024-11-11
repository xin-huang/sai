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
        x=0.5,
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

    assert df["Q95"].iloc[0] == 0.9, "Unexpected value in 'Q95' column"


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
        x=0.5,
        y=[("=", 1)],
        output_file=example_data["output_file"],
        stat_type="Q95",
        num_workers=1,
    )

    # Define output files for U and Q outliers
    u_outliers_file = example_data["output_dir"] / "U_outliers.tsv"
    q_outliers_file = example_data["output_dir"] / "Q_outliers.tsv"

    # Run the outlier function for U
    outlier(
        score_file=str(example_data["output_file"]),
        output=str(u_outliers_file),
        quantile=0.95,
    )

    # Run the outlier function for Q
    outlier(
        score_file=str(example_data["output_file"]),
        output=str(q_outliers_file),
        quantile=0.95,
    )

    # Check if the outlier files are created
    assert u_outliers_file.exists()
    assert q_outliers_file.exists()

    # Optionally, read and check contents of the outlier files
    u_outliers_df = pd.read_csv(u_outliers_file, sep="\t")
    q_outliers_df = pd.read_csv(q_outliers_file, sep="\t")

    assert u_outliers_df.empty, "U outliers file is unexpectedly not empty"
    assert q_outliers_df.empty, "Q outliers file is unexpectedly not empty"


@pytest.fixture
def sample_outlier_file():
    """Create a temporary outlier file to simulate input data."""
    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".tsv") as tmpfile:
        tmpfile.write("A\tB\tC\tD\tE\tF\tG\tH\tI\n")  # Example column names
        tmpfile.write("1\t2\t3\t4\t5\t6\t7\t0.1\t0.2\n")  # Example data
        tmpfile.write("2\t3\t4\t5\t6\t7\t8\t0.3\t0.4\n")
        tmpfile.write("3\t4\t5\t6\t7\t8\t9\t0.5\t0.6\n")
        return tmpfile.name


def test_plot(sample_outlier_file):
    """Test if the plot function correctly generates an output file."""
    with tempfile.NamedTemporaryFile(delete=False, suffix=".png") as tmp_output:
        output_file = tmp_output.name

    # Call the plot function
    plot(
        outlier_file=sample_outlier_file,
        output=output_file,
        xlabel="Q Values",
        ylabel="U Values",
        title="Test Plot",
    )

    # Check if the output file is created
    assert os.path.exists(output_file), "Plot output file was not created."

    # Clean up temporary files
    os.remove(sample_outlier_file)
    os.remove(output_file)
