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


import allel
import pytest
import numpy as np
from unittest.mock import mock_open, patch
from sai.utils import ChromosomeData
from sai.utils import check_anc_allele
from sai.utils import filter_fixed_variants
from sai.utils import filter_geno_data
from sai.utils import flip_snps
from sai.utils import get_ref_alt_allele
from sai.utils import parse_ind_file
from sai.utils import read_anc_allele
from sai.utils import read_data
from sai.utils import read_geno_data


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


# Test setup for read_geno_data
@pytest.fixture
def mock_vcf_data():
    return {
        "calldata/GT": np.array(
            [[[0, 1], [1, 1]], [[1, 0], [1, -1]], [[0, 0], [0, 0]]]
        ),
        "variants/CHROM": np.array(["chr1", "chr1", "chr2"]),
        "variants/POS": np.array([100, 200, 300]),
        "variants/REF": np.array(["A", "T", "G"]),
        "variants/ALT": np.array(["C", "A", "T"]),
        "samples": ["sample1", "sample2"],
    }


@patch("sai.utils.allel.read_vcf")
def test_read_geno_data(mock_read_vcf, mock_vcf_data):
    # Mock the return value of read_vcf
    mock_read_vcf.return_value = mock_vcf_data

    # Simulate parsed sample categories
    ind_samples = {"Category1": ["sample1"], "Category2": ["sample2"]}

    # Call the function with the mock data
    result = read_geno_data(
        vcf="mock.vcf",
        ind_samples=ind_samples,
        anc_allele_file=None,
        filter_missing=False,
    )

    # Verify the data structure for each chromosome
    assert "chr1" in result
    assert "chr2" in result

    # Check data shape and values for chr1 using ChromosomeData attributes
    chr1_data = result["chr1"]
    assert chr1_data.POS.tolist() == [100, 200]
    assert chr1_data.REF.tolist() == ["A", "T"]
    assert chr1_data.ALT.tolist() == ["C", "A"]
    assert chr1_data.GT.shape == (2, 2, 2)

    # Check data shape and values for chr2 using ChromosomeData attributes
    chr2_data = result["chr2"]
    assert chr2_data.POS.tolist() == [300]
    assert chr2_data.REF.tolist() == ["G"]
    assert chr2_data.ALT.tolist() == ["T"]
    assert chr2_data.GT.shape == (1, 2, 2)

    # Verify all samples are included
    assert result["samples"] == ["sample1", "sample2"]


# Test from files
@pytest.fixture
def data():
    pytest.ref_ind_list = "./tests/data/test.ref.ind.list"
    pytest.tgt_ind_list = "./tests/data/test.tgt.ind.list"
    pytest.vcf = "./tests/data/test.score.data.vcf"
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

    with pytest.raises(Exception) as e_info:
        emp_ind = parse_ind_file(pytest.emp_ind_list)


def test_read_geno_data_from_file(data):
    ref_ind = parse_ind_file(pytest.ref_ind_list)
    d = read_geno_data(pytest.vcf, ref_ind, None, filter_missing=False)

    vcf = allel.read_vcf(pytest.vcf, alt_number=1, samples=ref_ind["ref1"])

    assert np.array_equal(ref_ind["ref1"], vcf["samples"])
    assert np.array_equal(d["21"].POS, vcf["variants/POS"])
    assert np.array_equal(d["21"].REF, vcf["variants/REF"])
    assert np.array_equal(d["21"].ALT, vcf["variants/ALT"])
    assert np.array_equal(d["21"].GT, vcf["calldata/GT"])


def test_read_data_from_file(data):
    ref_data, ref_samples, tgt_data, tgt_samples, src_data, src_samples = read_data(
        vcf_file=pytest.vcf,
        ref_ind_file=pytest.ref_ind_list,
        tgt_ind_file=pytest.tgt_ind_list,
        src_ind_file=None,
        anc_allele_file=None,
        filter_ref=False,
        filter_tgt=False,
        filter_src=False,
    )

    rs = parse_ind_file(pytest.ref_ind_list)
    ts = parse_ind_file(pytest.tgt_ind_list)

    assert np.array_equal(rs, ref_samples)
    assert np.array_equal(ts, tgt_samples)

    ref_vcf = allel.read_vcf(pytest.vcf, alt_number=1, samples=rs["ref1"])
    tgt_vcf = allel.read_vcf(pytest.vcf, alt_number=1, samples=ts["tgt2"])

    assert np.array_equal(rs["ref1"], ref_vcf["samples"])
    assert np.array_equal(ts["tgt2"], tgt_vcf["samples"])
    assert np.array_equal(ref_data["ref1"]["21"].POS, ref_vcf["variants/POS"])
    assert np.array_equal(ref_data["ref1"]["21"].REF, ref_vcf["variants/REF"])
    assert np.array_equal(ref_data["ref1"]["21"].ALT, ref_vcf["variants/ALT"])
    assert np.array_equal(
        ref_data["ref1"]["21"].GT, ref_vcf["calldata/GT"].reshape(19, 4)
    )
    assert np.array_equal(tgt_data["tgt2"]["21"].POS, tgt_vcf["variants/POS"])
    assert np.array_equal(tgt_data["tgt2"]["21"].REF, tgt_vcf["variants/REF"])
    assert np.array_equal(tgt_data["tgt2"]["21"].ALT, tgt_vcf["variants/ALT"])
    assert np.array_equal(
        tgt_data["tgt2"]["21"].GT, tgt_vcf["calldata/GT"].reshape(19, 4)
    )


def test_read_anc_allele(data):
    anc_allele = read_anc_allele(pytest.anc_allele)

    exp_anc_allele = {"21": {2309: "G", 7879: "A", 11484: "-", 48989: "C"}}

    assert anc_allele == exp_anc_allele

    with pytest.raises(Exception) as e_info:
        anc_allele = read_anc_allele(pytest.emp_anc_allele)


def test_get_ref_alt_allele(data):
    tgt_ind = parse_ind_file(pytest.tgt_ind_list)
    tgt_vcf = allel.read_vcf(pytest.vcf, alt_number=1, samples=tgt_ind["tgt1"])

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
    ref_data, ref_samples, tgt_data, tgt_samples, src_data, src_samples = read_data(
        vcf_file=pytest.vcf,
        ref_ind_file=pytest.ref_ind_list,
        tgt_ind_file=pytest.tgt_ind_list,
        src_ind_file=None,
        anc_allele_file=pytest.anc_allele,
        filter_ref=False,
        filter_tgt=False,
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

    assert np.array_equal(ref_data["ref1"]["21"].GT, exp_ref_gt.reshape(3, 4))
    assert np.array_equal(tgt_data["tgt1"]["21"].GT, exp_tgt_gt1.reshape(3, 4))
    assert np.array_equal(tgt_data["tgt2"]["21"].GT, exp_tgt_gt2.reshape(3, 4))
    assert np.array_equal(tgt_data["tgt1"]["21"].POS, exp_tgt_pos)
    assert np.array_equal(tgt_data["tgt2"]["21"].POS, exp_tgt_pos)


# Test data setup for filter_fixed_variants
@pytest.fixture
def sample_data():
    # Sample ChromosomeData with mixed fixed and non-fixed variants
    return {
        "pop1": {
            "chr1": ChromosomeData(
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
            "chr2": ChromosomeData(
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
    }


@pytest.fixture
def sample_info():
    # Sample information for two individuals in 'pop1'
    return {"pop1": ["sample1", "sample2"]}


def test_filter_fixed_variants(sample_data, sample_info):
    # Apply the filter_fixed_variants function
    filtered_data = filter_fixed_variants(sample_data, sample_info)

    # Verify that fixed variants are removed
    assert "pop1" in filtered_data
    assert "chr1" in filtered_data["pop1"]
    assert "chr2" in filtered_data["pop1"]

    # Check that only non-fixed variants are retained for chr1
    chr1_data = filtered_data["pop1"]["chr1"]
    assert chr1_data.POS.tolist() == [300, 400]  # Positions with mixed genotypes
    assert chr1_data.REF.tolist() == ["T", "C"]
    assert chr1_data.ALT.tolist() == ["G", "T"]
    assert chr1_data.GT.shape == (2, 2, 2)

    # Verify that all variants in chr2 are filtered out, as they are fixed
    chr2_data = filtered_data["pop1"]["chr2"]
    assert chr2_data.POS.size == 0
    assert chr2_data.REF.size == 0
    assert chr2_data.ALT.size == 0
    assert chr2_data.GT.size == 0


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
