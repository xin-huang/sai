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
import numpy as np
from sai.stats import calc_freq
from sai.stats import compute_matching_loci
from sai.stats import calc_four_pops_freq
from sai.stats import calc_pattern_sum


def test_phased_data():
    # Phased data, ploidy = 1
    gts = np.array([[1, 0, 0, 1], [0, 0, 0, 0], [1, 1, 1, 1]])
    expected_frequency = np.array([0.5, 0.0, 1.0])
    result = calc_freq(gts, ploidy=1)
    np.testing.assert_array_almost_equal(
        result, expected_frequency, decimal=6, err_msg="Phased data test failed."
    )


def test_unphased_diploid_data():
    # Unphased data, ploidy = 2 (diploid)
    gts = np.array([[1, 1], [0, 0], [2, 2]])
    expected_frequency = np.array([0.5, 0.0, 1.0])
    result = calc_freq(gts, ploidy=2)
    np.testing.assert_array_almost_equal(
        result,
        expected_frequency,
        decimal=6,
        err_msg="Unphased diploid data test failed.",
    )


def test_unphased_triploid_data():
    # Unphased data, ploidy = 3 (triploid)
    gts = np.array([[1, 2, 3], [0, 0, 0], [3, 3, 3]])
    expected_frequency = np.array([0.6667, 0.0, 1.0])
    result = calc_freq(gts, ploidy=3)
    np.testing.assert_array_almost_equal(
        result,
        expected_frequency,
        decimal=4,
        err_msg="Unphased triploid data test failed.",
    )


def test_unphased_tetraploid_data():
    # Unphased data, ploidy = 4 (tetraploid)
    gts = np.array([[2, 2, 2, 2], [1, 3, 0, 4], [0, 0, 0, 0]])
    expected_frequency = np.array([0.5, 0.5, 0.0])
    result = calc_freq(gts, ploidy=4)
    np.testing.assert_array_almost_equal(
        result,
        expected_frequency,
        decimal=6,
        err_msg="Unphased tetraploid data test failed.",
    )


def test_invalid_ploidy():
    gts = np.array([[1, 2, 3], [0, 0, 0], [3, 3, 3]])

    with pytest.raises(ValueError):
        calc_freq(gts, ploidy=None)

    with pytest.raises(ValueError):
        calc_freq(gts, ploidy=9.9)

    with pytest.raises(ValueError):
        calc_freq(gts, ploidy=-100)


def test_compute_matching_loci():
    # Sample genotype data
    ref_gts = np.array([[0, 1, 0], [1, 1, 0], [0, 0, 1]])
    tgt_gts = np.array([[1, 1, 0], [0, 1, 1], [1, 1, 1]])
    src_gts_list = [
        np.array([[0, 0, 1], [1, 1, 0], [0, 1, 1]]),  # src1
        np.array([[1, 1, 0], [1, 0, 0], [1, 1, 0]]),  # src2
    ]

    # Define parameters with all possible conditions
    conditions = [("=", 0.5), ("<", 0.4), (">", 0.3), ("<=", 0.6), (">=", 0.2)]
    ploidy = [2, 2, 2]
    anc_allele_available = False

    for y_condition in conditions:
        y_list = [y_condition, y_condition]  # Apply the same condition to both sources

        # Call the function
        ref_freq, tgt_freq, condition = compute_matching_loci(
            ref_gts,
            tgt_gts,
            src_gts_list,
            0.5,
            y_list,
            ploidy,
            anc_allele_available,
        )

        # Assertions to verify the outputs
        assert ref_freq.shape == (3,)
        assert tgt_freq.shape == (3,)
        assert condition.shape == (3,)
        assert np.all((ref_freq >= 0) & (ref_freq <= 1))
        assert np.all((tgt_freq >= 0) & (tgt_freq <= 1))
        assert np.all(
            np.logical_or(condition == True, condition == False)
        )  # Ensure condition is boolean

    # Test invalid w values
    with pytest.raises(
        ValueError, match=r"Parameters w must be within the range \[0, 1\]."
    ):
        compute_matching_loci(
            ref_gts,
            tgt_gts,
            src_gts_list,
            -0.1,
            y_list,
            ploidy,
            anc_allele_available,
        )
    with pytest.raises(
        ValueError, match=r"Parameters w must be within the range \[0, 1\]."
    ):
        compute_matching_loci(
            ref_gts,
            tgt_gts,
            src_gts_list,
            1.1,
            y_list,
            ploidy,
            anc_allele_available,
        )

    # Test invalid y values
    with pytest.raises(ValueError, match="Invalid value in y_list"):
        compute_matching_loci(
            ref_gts,
            tgt_gts,
            src_gts_list,
            0.5,
            [("=", -0.1)],
            ploidy,
            anc_allele_available,
        )
    with pytest.raises(ValueError, match="Invalid value in y_list"):
        compute_matching_loci(
            ref_gts,
            tgt_gts,
            src_gts_list,
            0.5,
            [("=", 1.1)],
            ploidy,
            anc_allele_available,
        )

    # Test invalid operators
    with pytest.raises(ValueError, match="Invalid operator in y_list"):
        compute_matching_loci(
            ref_gts,
            tgt_gts,
            src_gts_list,
            0.5,
            [("invalid", 0.5)],
            ploidy,
            anc_allele_available,
        )

    # Test mismatched src_gts_list and y_list lengths
    with pytest.raises(
        ValueError, match="The length of src_gts_list and y_list must match"
    ):
        compute_matching_loci(
            ref_gts,
            tgt_gts,
            src_gts_list,
            0.5,
            [("=", 0.5)],
            ploidy,
            anc_allele_available,
        )


def test_calc_four_pops_freq_basic():
    ref_gts = np.array([[0, 1], [1, 1]])  # freq = [0.5, 1.0]
    tgt_gts = np.array([[1, 0], [0, 0]])  # freq = [0.5, 0.0]
    src_gts = np.array([[1, 1], [1, 0]])  # freq = [1.0, 0.5]
    out_gts = np.array([[0, 0], [0, 1]])  # freq = [0.0, 0.5]

    ref, tgt, src, out = calc_four_pops_freq(
        ref_gts,
        tgt_gts,
        src_gts,
        out_gts,
    )

    np.testing.assert_array_almost_equal(ref, np.array([0.5, 1.0]))
    np.testing.assert_array_almost_equal(tgt, np.array([0.5, 0.0]))
    np.testing.assert_array_almost_equal(src, np.array([1.0, 0.5]))
    np.testing.assert_array_almost_equal(out, np.array([0.0, 0.5]))


def test_calc_four_pops_freq_no_outgroup():
    ref_gts = np.array([[0, 1]])
    tgt_gts = np.array([[1, 0]])
    src_gts = np.array([[1, 1]])

    ref, tgt, src, out = calc_four_pops_freq(ref_gts, tgt_gts, src_gts, out_gts=None)

    np.testing.assert_array_equal(ref, np.array([0.5]))
    np.testing.assert_array_equal(tgt, np.array([0.5]))
    np.testing.assert_array_equal(src, np.array([1.0]))
    np.testing.assert_array_equal(out, np.array([0.0]))  # default to 0s


def test_calc_four_pops_freq_diploid():
    ref_gts = np.array([[0, 2]])
    tgt_gts = np.array([[1, 1]])
    src_gts = np.array([[2, 0]])
    out_gts = np.array([[1, 1]])

    # ploidy=2 → total alleles = 2 * n_samples
    # freq = sum / (2 * N)

    ref, tgt, src, out = calc_four_pops_freq(
        ref_gts=ref_gts,
        tgt_gts=tgt_gts,
        src_gts=src_gts,
        out_gts=out_gts,
        ref_ploidy=2,
        tgt_ploidy=2,
        src_ploidy=2,
        out_ploidy=2,
    )

    np.testing.assert_array_equal(ref, np.array([0.5]))  # (0+2)/4
    np.testing.assert_array_equal(tgt, np.array([0.5]))  # (1+1)/4
    np.testing.assert_array_equal(src, np.array([0.5]))  # (2+0)/4
    np.testing.assert_array_equal(out, np.array([0.5]))  # (1+1)/4


def test_calc_four_pops_freq_mixed_ploidy():
    ref_gts = np.array([[0, 2]])
    tgt_gts = np.array([[1, 1]])
    src_gts = np.array([[2, 0]])
    out_gts = np.array([[1, 1]])

    ref, tgt, src, out = calc_four_pops_freq(
        ref_gts=ref_gts,
        tgt_gts=tgt_gts,
        src_gts=src_gts,
        out_gts=out_gts,
        ref_ploidy=2,
        tgt_ploidy=1,
        src_ploidy=4,
        out_ploidy=4,
    )

    np.testing.assert_array_equal(ref, np.array([0.5]))  # (0+2)/4
    np.testing.assert_array_equal(tgt, np.array([1]))  # (1+1)/2
    np.testing.assert_array_equal(src, np.array([0.25]))  # (2+0)/8
    np.testing.assert_array_equal(out, np.array([0.25]))  # (1+1)/8


def test_calc_pattern_sum_abba():
    ref = np.array([0.1, 0.8])
    tgt = np.array([0.9, 0.2])
    src = np.array([0.5, 0.5])
    out = np.array([0.0, 1.0])

    # pattern: 'abba'
    # site 0: (1-0.1)*0.9*0.5*(1-0.0) = 0.9*0.9*0.5*1 = 0.405
    # site 1: (1-0.8)*0.2*0.5*(1-1.0) = 0.2*0.2*0.5*0 = 0.0
    # sum = 0.405 + 0.0 = 0.405

    result = calc_pattern_sum(ref, tgt, src, out, "abba")
    assert np.isclose(result, 0.405)


def test_calc_pattern_sum_baba():
    ref = np.array([0.1, 0.8])
    tgt = np.array([0.9, 0.2])
    src = np.array([0.5, 0.5])
    out = np.array([0.0, 1.0])

    # pattern: 'baba'
    # site 0: 0.1*(1-0.9)*0.5*(1-0.0) = 0.1*0.1*0.5*1 = 0.005
    # site 1: 0.8*(1-0.2)*0.5*0 = 0.8*0.8*0.5*0 = 0
    # sum = 0.005

    result = calc_pattern_sum(ref, tgt, src, out, "baba")
    assert np.isclose(result, 0.005)


def test_calc_pattern_sum_baaa():
    ref = np.array([0.1, 0.8])
    tgt = np.array([0.9, 0.2])
    src = np.array([0.5, 0.5])
    out = np.array([0.0, 1.0])

    # pattern: 'baaa'
    # site 0: 0.1*(1-0.9)*(1-0.5)*(1-0.0) = 0.1*0.1*0.5*1 = 0.005
    # site 1: 0.8*(1-0.2)*(1-0.5)*0      = 0.8*0.8*0.5*0 = 0
    # sum = 0.005

    result = calc_pattern_sum(ref, tgt, src, out, "baaa")
    assert np.isclose(result, 0.005)


def test_calc_pattern_sum_abaa():
    ref = np.array([0.1, 0.8])
    tgt = np.array([0.9, 0.2])
    src = np.array([0.5, 0.5])
    out = np.array([0.0, 1.0])

    # pattern: 'abaa'
    # site 0: (1-0.1)*0.9*(1-0.5)*(1-0.0) = 0.9*0.9*0.5*1 = 0.405
    # site 1: (1-0.8)*0.2*(1-0.5)*0       = 0.2*0.2*0.5*0 = 0
    # sum = 0.405

    result = calc_pattern_sum(ref, tgt, src, out, "abaa")
    assert np.isclose(result, 0.405)


def test_invalid_pattern_length():
    ref = tgt = src = out = np.array([0.1, 0.2])
    with pytest.raises(ValueError, match="four-character"):
        _ = calc_pattern_sum(ref, tgt, src, out, "ab")


def test_invalid_pattern_char():
    ref = tgt = src = out = np.array([0.1, 0.2])
    with pytest.raises(ValueError, match="Invalid character"):
        _ = calc_pattern_sum(ref, tgt, src, out, "abxa")
