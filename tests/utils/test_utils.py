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


import allel
import pytest
import numpy as np
import pandas as pd
from unittest.mock import mock_open, patch
from sai.configs import PloidyConfig
from sai.utils import ChromosomeData
from sai.utils import filter_fixed_variants
from sai.utils import filter_geno_data
from sai.utils import flip_snps
from sai.utils import get_ref_alt_allele
from sai.utils import parse_ind_file
from sai.utils import read_anc_allele
from sai.utils import read_data
from sai.utils import read_geno_data
from sai.utils import split_genome
from sai.utils import natsorted_df


def test_valid_file():
    # Mock the content of a valid file with categories
    mock_data = "Category1 Sample1\nCategory1 Sample2\nCategory2 Sample3\n"
    with patch("builtins.open", mock_open(read_data=mock_data)):
        samples = parse_ind_file("mock_file.txt")
        assert samples == {
            "Category1": ["Sample1", "Sample2"],
            "Category2": ["Sample3"],
        }


def test_empty_file():
    # Mock an empty file
    mock_data = ""
    with patch("builtins.open", mock_open(read_data=mock_data)):
        with pytest.raises(ValueError) as excinfo:
            parse_ind_file("mock_file.txt")
        assert (
            str(excinfo.value)
            == "No samples found in mock_file.txt. Please check your data."
        )


def test_file_not_found():
    # Ensure FileNotFoundError is raised when file does not exist
    with pytest.raises(FileNotFoundError):
        parse_ind_file("non_existent_file.txt")


def test_ignores_empty_lines():
    # Mock a file with empty lines and valid lines
    mock_data = "Category1 Sample1\n\nCategory1 Sample2\n   \nCategory2 Sample3\n"
    with patch("builtins.open", mock_open(read_data=mock_data)):
        samples = parse_ind_file("mock_file.txt")
        assert samples == {
            "Category1": ["Sample1", "Sample2"],
            "Category2": ["Sample3"],
        }


# Test data setup for filter_geno_data
@pytest.fixture
def sample_genotype_data():
    return ChromosomeData(
        POS=np.array([100, 200, 300, 400, 500]),
        REF=np.array(["A", "T", "G", "C", "A"]),
        ALT=np.array(["C", "A", "T", "G", "T"]),
        GT=np.array(
            [
                [[0, 1], [1, 1]],
                [[1, 0], [1, -1]],
                [[0, 0], [0, 0]],
                [[1, 1], [1, 1]],
                [[0, 1], [0, 0]],
            ]
        ),
    )


def test_filter_geno_data(sample_genotype_data):
    # Example of filtering out the third row
    index = np.array([True, False, True, True, False])
    filtered = filter_geno_data(sample_genotype_data, index)

    # Assertions updated to check ChromosomeData attributes
    assert filtered.POS.tolist() == [100, 300, 400]
    assert filtered.REF.tolist() == ["A", "G", "C"]
    assert filtered.ALT.tolist() == ["C", "T", "G"]
    assert filtered.GT.shape == (3, 2, 2)


# Test from files
@pytest.fixture
def data():
    pytest.ref_ind_list = "./tests/data/test.ref.ind.list"
    pytest.tgt_ind_list = "./tests/data/test.tgt.ind.list"
    pytest.vcf = "./tests/data/test.data.vcf"
    pytest.anc_allele = "./tests/data/test.anc.allele.bed"


def test_parse_ind_file_from_files(data):
    ref_ind = parse_ind_file(pytest.ref_ind_list)
    tgt_ind = parse_ind_file(pytest.tgt_ind_list)

    exp_ref_ind = {
        "ref1": ["ind5", "ind6"],
    }
    exp_tgt_ind = {
        "tgt1": ["ind1", "ind2"],
        "tgt2": ["ind3", "ind4"],
    }

    assert ref_ind == exp_ref_ind
    assert tgt_ind == exp_tgt_ind


def test_read_geno_data_from_file(data):
    ref_ind = parse_ind_file(pytest.ref_ind_list)
    d = read_geno_data(
        vcf=pytest.vcf,
        ind_samples=ref_ind,
        chr_name="21",
        anc_allele_file=None,
        filter_missing=False,
    )

    vcf = allel.read_vcf(pytest.vcf, alt_number=1, samples=ref_ind["ref1"], region="21")

    assert np.array_equal(ref_ind["ref1"], vcf["samples"])
    assert np.array_equal(d["ref1"].POS, vcf["variants/POS"])
    assert np.array_equal(d["ref1"].REF, vcf["variants/REF"])
    assert np.array_equal(d["ref1"].ALT, vcf["variants/ALT"])
    assert np.array_equal(d["ref1"].GT, vcf["calldata/GT"])


def test_read_data_from_file(data):
    ploidy_config = PloidyConfig(
        {
            "ref": {"ref1": 2},
            "tgt": {"tgt1": 2, "tgt2": 2},
            "src": {"src1": 2, "src2": 2},
        }
    )

    results = read_data(
        vcf_file=pytest.vcf,
        chr_name="21",
        ref_ind_file=pytest.ref_ind_list,
        tgt_ind_file=pytest.tgt_ind_list,
        src_ind_file=None,
        out_ind_file=None,
        anc_allele_file=None,
        filter_ref=False,
        filter_tgt=False,
        filter_src=False,
        filter_out=False,
        ploidy_config=ploidy_config,
    )

    rs = parse_ind_file(pytest.ref_ind_list)
    ts = parse_ind_file(pytest.tgt_ind_list)

    assert np.array_equal(rs, results["ref"][1])
    assert np.array_equal(ts, results["tgt"][1])

    ref_vcf = allel.read_vcf(pytest.vcf, alt_number=1, samples=rs["ref1"], region="21")
    tgt_vcf = allel.read_vcf(pytest.vcf, alt_number=1, samples=ts["tgt2"], region="21")

    assert np.array_equal(rs["ref1"], ref_vcf["samples"])
    assert np.array_equal(ts["tgt2"], tgt_vcf["samples"])
    assert np.array_equal(results["ref"][0]["ref1"].POS, ref_vcf["variants/POS"])
    assert np.array_equal(results["ref"][0]["ref1"].REF, ref_vcf["variants/REF"])
    assert np.array_equal(results["ref"][0]["ref1"].ALT, ref_vcf["variants/ALT"])
    assert np.array_equal(
        results["ref"][0]["ref1"].GT, ref_vcf["calldata/GT"].reshape(19, 4)
    )
    assert np.array_equal(results["tgt"][0]["tgt2"].POS, tgt_vcf["variants/POS"])
    assert np.array_equal(results["tgt"][0]["tgt2"].REF, tgt_vcf["variants/REF"])
    assert np.array_equal(results["tgt"][0]["tgt2"].ALT, tgt_vcf["variants/ALT"])
    assert np.array_equal(
        results["tgt"][0]["tgt2"].GT, tgt_vcf["calldata/GT"].reshape(19, 4)
    )


def test_read_anc_allele(data):
    anc_allele = read_anc_allele(pytest.anc_allele, "21")

    exp_anc_allele = {"21": {2309: "G", 7879: "A", 11484: "-", 48989: "C"}}

    assert anc_allele == exp_anc_allele


def test_get_ref_alt_allele(data):
    tgt_ind = parse_ind_file(pytest.tgt_ind_list)
    tgt_vcf = allel.read_vcf(
        pytest.vcf, alt_number=1, samples=tgt_ind["tgt1"], region="21"
    )

    ref_allele, alt_allele = get_ref_alt_allele(
        tgt_vcf["variants/REF"], tgt_vcf["variants/ALT"], tgt_vcf["variants/POS"]
    )

    exp_ref_allele = {
        2309: "G",
        7879: "C",
        11484: "A",
        16249: "A",
        17324: "G",
        19064: "G",
        19124: "G",
        23559: "G",
        25354: "G",
        26654: "G",
        29724: "G",
        30769: "C",
        31319: "C",
        37199: "C",
        38009: "C",
        39444: "C",
        40809: "C",
        45079: "C",
        48989: "C",
    }
    exp_alt_allele = {
        2309: "A",
        7879: "A",
        11484: "C",
        16249: "C",
        17324: "T",
        19064: "T",
        19124: "A",
        23559: "A",
        25354: "T",
        26654: "C",
        29724: "A",
        30769: "T",
        31319: "T",
        37199: "T",
        38009: "T",
        39444: "T",
        40809: "T",
        45079: "T",
        48989: "T",
    }

    assert ref_allele == exp_ref_allele
    assert alt_allele == exp_alt_allele


def test_check_anc_allele(data):
    ploidy_config = PloidyConfig(
        {
            "ref": {"ref1": 2},
            "tgt": {"tgt1": 2, "tgt2": 2},
            "src": {"src1": 2, "src2": 2},
        }
    )

    data = read_data(
        vcf_file=pytest.vcf,
        chr_name="21",
        ref_ind_file=pytest.ref_ind_list,
        tgt_ind_file=pytest.tgt_ind_list,
        src_ind_file=None,
        out_ind_file=None,
        anc_allele_file=pytest.anc_allele,
        filter_ref=False,
        filter_tgt=False,
        ploidy_config=ploidy_config,
    )

    exp_ref_gt = allel.GenotypeArray(
        [
            [[0, 0], [0, 0]],
            [[1, 1], [1, 1]],
            [[0, 0], [0, 0]],
        ],
    )
    exp_tgt_gt1 = allel.GenotypeArray(
        [
            [[1, 0], [0, 0]],
            [[1, 1], [1, 0]],
            [[0, 0], [0, 1]],
        ],
    )
    exp_tgt_gt2 = allel.GenotypeArray(
        [
            [[0, 0], [0, 0]],
            [[1, 1], [1, 1]],
            [[0, 0], [0, 0]],
        ],
    )
    exp_tgt_pos = [2309, 7879, 48989]

    assert np.array_equal(data["ref"][0]["ref1"].GT, exp_ref_gt.reshape(3, 4))
    assert np.array_equal(data["tgt"][0]["tgt1"].GT, exp_tgt_gt1.reshape(3, 4))
    assert np.array_equal(data["tgt"][0]["tgt2"].GT, exp_tgt_gt2.reshape(3, 4))
    assert np.array_equal(data["tgt"][0]["tgt1"].POS, exp_tgt_pos)
    assert np.array_equal(data["tgt"][0]["tgt2"].POS, exp_tgt_pos)


# Test data setup for filter_fixed_variants
@pytest.fixture
def sample_data():
    # Sample ChromosomeData with mixed fixed and non-fixed variants
    return {
        "pop1": ChromosomeData(
            POS=np.array([100, 200, 300, 400]),
            REF=np.array(["A", "G", "T", "C"]),
            ALT=np.array(["C", "A", "G", "T"]),
            GT=allel.GenotypeArray(
                [
                    [[0, 0], [0, 0]],  # Fixed ref (AA, AA)
                    [[1, 1], [1, 1]],  # Fixed alt (CC, CC)
                    [[0, 1], [1, 1]],  # Mixed (AG, GG)
                    [[0, 1], [0, 0]],  # Mixed (AC, AA)
                ]
            ),
        ),
        "pop2": ChromosomeData(
            POS=np.array([150, 250]),
            REF=np.array(["T", "A"]),
            ALT=np.array(["G", "C"]),
            GT=allel.GenotypeArray(
                [
                    [[0, 0], [0, 0]],  # Fixed ref (TT, TT)
                    [[1, 1], [1, 1]],  # Fixed alt (CC, CC)
                ]
            ),
        ),
    }


@pytest.fixture
def sample_info():
    # Sample information for two individuals in 'pop1'
    return {
        "pop1": ["sample1", "sample2"],
        "pop2": ["sample3", "sample4"],
    }


def test_filter_fixed_variants(sample_data, sample_info):
    # Apply the filter_fixed_variants function
    filtered_data = filter_fixed_variants(sample_data, sample_info)

    # Verify that fixed variants are removed
    assert "pop1" in filtered_data
    assert "pop2" in filtered_data

    # Check that only non-fixed variants are retained for chr1
    pop1_data = filtered_data["pop1"]
    assert pop1_data.POS.tolist() == [300, 400]  # Positions with mixed genotypes
    assert pop1_data.REF.tolist() == ["T", "C"]
    assert pop1_data.ALT.tolist() == ["G", "T"]
    assert pop1_data.GT.shape == (2, 2, 2)

    # Verify that all variants in chr2 are filtered out, as they are fixed
    pop2_data = filtered_data["pop2"]
    assert pop2_data.POS.size == 0
    assert pop2_data.REF.size == 0
    assert pop2_data.ALT.size == 0
    assert pop2_data.GT.size == 0


# Test data setup for flip_snps
@pytest.fixture
def sample_chromosome_data():
    # Sample ChromosomeData with genotypes to test flipping
    return ChromosomeData(
        POS=np.array([100, 200, 300, 400]),
        REF=np.array(["A", "G", "T", "C"]),
        ALT=np.array(["C", "A", "G", "T"]),
        GT=allel.GenotypeArray(
            [
                [[0, 1], [1, 1]],  # Mixed genotype, should be flipped
                [[1, 1], [1, 1]],  # ALT fixed, should be flipped
                [[0, 0], [0, 1]],  # Mixed genotype, should remain unchanged
                [[1, 0], [0, 0]],  # Mixed genotype, should be flipped
            ]
        ),
    )


def test_flip_snps(sample_chromosome_data):
    # Define SNP positions to be flipped
    flipped_snps = [100, 200, 400]

    # Apply the flip_snps function
    flip_snps(sample_chromosome_data, flipped_snps)

    # Check flipped genotypes
    # Position 100: original [[0, 1], [1, 1]] -> flipped [[1, 0], [0, 0]]
    assert sample_chromosome_data.GT[0].tolist() == [[1, 0], [0, 0]]
    # Position 200: original [[1, 1], [1, 1]] -> flipped [[0, 0], [0, 0]]
    assert sample_chromosome_data.GT[1].tolist() == [[0, 0], [0, 0]]
    # Position 300: original [[0, 0], [0, 1]] -> should remain unchanged
    assert sample_chromosome_data.GT[2].tolist() == [[0, 0], [0, 1]]
    # Position 400: original [[1, 0], [0, 0]] -> flipped [[0, 1], [1, 1]]
    assert sample_chromosome_data.GT[3].tolist() == [[0, 1], [1, 1]]


# Sample test function for split_genome
def test_split_genome():
    # Test case 1: Basic case with regular windows
    pos = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    window_size = 30
    step_size = 20
    result = split_genome(pos, window_size, step_size)
    expected = [(1, 30), (21, 50), (41, 70), (61, 90), (81, 110)]
    assert result == expected, f"Expected {expected}, but got {result}"

    # Test case 2: Step size is larger than window size
    pos = np.array([0, 10, 20, 30, 40, 50])
    window_size = 20
    step_size = 25

    with pytest.raises(
        ValueError, match="`step_size` cannot be greater than `window_size`"
    ):
        split_genome(pos, window_size, step_size)

    # Test case 3: Handle empty `pos` array
    pos = np.array([])
    window_size = 30
    step_size = 10
    with pytest.raises(ValueError, match="`pos` array must not be empty"):
        split_genome(pos, window_size, step_size)


def test_natsorted_df_correct_order():
    df = pd.DataFrame(
        {
            "Chrom": ["1", "10", "2", "X", "1"],
            "Start": [300, 50, 150, 10, 100],
            "End": [400, 100, 200, 50, 200],
        }
    )

    sorted_df = natsorted_df(df)

    expected_df = pd.DataFrame(
        {
            "Chrom": ["1", "1", "2", "10", "X"],
            "Start": [100, 300, 150, 50, 10],
            "End": [200, 400, 200, 100, 50],
        }
    ).reset_index(drop=True)

    pd.testing.assert_frame_equal(sorted_df, expected_df)


def test_natsorted_df_missing_columns():
    df_missing = pd.DataFrame(
        {
            "Chrom": ["1", "2", "X"],
            "Start": [100, 200, 300],
        }
    )

    with pytest.raises(ValueError, match="Missing required columns: End"):
        natsorted_df(df_missing)


def test_natsorted_df_empty_dataframe():
    df_empty = pd.DataFrame(columns=["Chrom", "Start", "End"])
    sorted_df = natsorted_df(df_empty)

    assert sorted_df.empty


def test_natsorted_df_single_row():
    df_single = pd.DataFrame({"Chrom": ["1"], "Start": [100], "End": [200]})

    sorted_df = natsorted_df(df_single)

    pd.testing.assert_frame_equal(sorted_df, df_single)


def test_natsorted_df_integer_start_end():
    df_mixed_types = pd.DataFrame(
        {
            "Chrom": ["1", "2", "X"],
            "Start": ["100", "200", "300"],
            "End": ["150", "250", "350"],
        }
    )

    sorted_df = natsorted_df(df_mixed_types)

    assert sorted_df["Start"].dtype == int
    assert sorted_df["End"].dtype == int
