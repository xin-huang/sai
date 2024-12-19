# Copyright 2024 Xin Huang and Florian R. Schmidt
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
from sai.stats.features import calc_u
from sai.stats.features import calc_q
from sai.stats.features import calc_freq
from sai.stats.features import calc_seq_div
from sai.stats.features import calc_rd
from sai.stats.features import calc_abba_baba


def test_calc_u_basic():
    # Test data
    ref_gts = np.array([[0, 0, 1], [0, 0, 0], [1, 1, 1]])
    tgt_gts = np.array([[1, 1, 1], [1, 0, 0], [0, 1, 0]])
    src_gts = np.array([[0, 0, 0], [1, 1, 1], [1, 0, 1]])

    # Parameters w, x, y values
    w, x, y = 0.5, 0.5, 0

    # Expected output
    expected_count = 1  # Only the first site meets the criteria

    # Run test
    result = calc_u(ref_gts, tgt_gts, [src_gts], w, x, [y])
    assert result == expected_count


def test_calc_u_no_match():
    # Test data with no matching sites
    ref_gts = np.array([[0, 1, 1], [1, 1, 1]])
    tgt_gts = np.array([[0, 0, 0], [1, 0, 1]])
    src_gts = np.array([[1, 1, 1], [1, 1, 1]])

    # Parameters for testing with no matching sites
    w, x, y = 0.3, 0.5, 0

    expected_count = 0  # No sites meet the criteria

    result = calc_u(ref_gts, tgt_gts, [src_gts], w, x, [y])
    assert result == expected_count


def test_calc_u_all_match():
    # Test data where all sites meet the criteria
    ref_gts = np.array([[0, 0, 0], [0, 0, 0]])
    tgt_gts = np.array([[1, 1, 1], [1, 1, 1]])
    src_gts = np.array([[0, 0, 0], [0, 0, 0]])

    # Parameters for testing all sites matching
    w, x, y = 0.5, 0.5, 0

    expected_count = 2  # All sites meet the criteria

    result = calc_u(ref_gts, tgt_gts, [src_gts], w, x, [y])
    assert result == expected_count


def test_calc_q_basic():
    # Test data
    ref_gts = np.array([[0, 0, 1], [0, 0, 0], [1, 1, 1]])
    tgt_gts = np.array([[0, 1, 1], [0, 0, 1], [1, 1, 1]])
    src_gts = np.array([[1, 1, 1], [0, 1, 1], [1, 1, 1]])
    w, y, quantile = 0.5, 1.0, 0.95

    # Expected output
    expected_result = 0.66667  # Only the first site meets the criteria

    # Run test
    result = calc_q(ref_gts, tgt_gts, [src_gts], w, [y], quantile)
    assert np.isclose(
        result, expected_result
    ), f"Expected {expected_result}, got {result}"


def test_calc_q_no_match():
    # Test data with no matching loci
    ref_gts = np.array([[0, 0, 1], [1, 1, 1]])
    tgt_gts = np.array([[0, 1, 1], [1, 1, 1]])
    src_gts = np.array([[1, 1, 1], [1, 1, 1]])
    w, y, quantile = (
        0.3,
        0.0,
        0.95,
    )  # No tgt_gts frequencies < w and no src_gts frequencies == y

    # Expected output
    expected_result = np.nan  # No loci meet criteria, should return NaN

    # Run test
    result = calc_q(ref_gts, tgt_gts, [src_gts], w, [y], quantile)
    assert np.isnan(result), f"Expected NaN, got {result}"


def test_calc_q_different_quantile():
    # Test data
    ref_gts = np.array([[0, 0, 1], [1, 0, 0], [0, 0, 1]])
    tgt_gts = np.array([[0, 1, 1], [1, 1, 1], [1, 1, 1]])
    src_gts = np.array([[0, 0, 0], [1, 1, 1], [1, 1, 1]])
    w, y, quantile = 0.5, 1.0, 0.5

    # Expected output
    expected_result = (
        1.0  # 50% quantile (median) of [1.0, 1.0, 1.0] in tgt_gts that meets conditions
    )

    # Run test
    result = calc_q(ref_gts, tgt_gts, [src_gts], w, [y], quantile)
    assert np.isclose(
        result, expected_result
    ), f"Expected {expected_result}, got {result}"


def test_calc_q_edge_case():
    # Edge case where only one site meets criteria
    ref_gts = np.array([[0, 0, 1], [0, 0, 0], [1, 1, 1]])
    tgt_gts = np.array([[0, 1, 1], [1, 1, 1], [0, 0, 0]])
    src_gts = np.array([[0, 0, 0], [1, 1, 1], [1, 1, 1]])
    w, y, quantile = 0.95, 1.0, 0.95

    # Expected output
    expected_result = 1.0  # Only one matching site in tgt_gts

    # Run test
    result = calc_q(ref_gts, tgt_gts, [src_gts], w, [y], quantile)
    assert np.isclose(
        result, expected_result
    ), f"Expected {expected_result}, got {result}"


def test_calc_q_invalid_quantile():
    # Test data for invalid quantile
    ref_gts = np.array([[0, 0, 1]])
    tgt_gts = np.array([[0, 1, 1]])
    src_gts = np.array([[1, 1, 1]])
    w, y, quantile = 0.5, 1.0, 1.5  # Invalid quantile (out of [0, 1] range)

    with pytest.raises(
        ValueError,
        match="Parameters w and quantile must be within the range \\[0, 1\\]",
    ):
        calc_q(ref_gts, tgt_gts, [src_gts], w, [y], quantile)


def test_calc_u_with_two_sources():
    # Test data
    ref_gts = np.array([[0, 0, 1], [0, 0, 0], [1, 1, 1]])
    tgt_gts = np.array([[0, 1, 1], [0, 0, 1], [1, 1, 1]])
    src_gts1 = np.array([[1, 1, 1], [0, 1, 1], [1, 1, 1]])
    src_gts2 = np.array([[1, 1, 1], [0, 1, 1], [1, 1, 1]])
    w, x, y_list = 0.5, 0.5, [1.0, 1.0]

    # Expected result: only loci where ref_freq < w, tgt_freq > x, src_freq1 == y1, and src_freq2 == y2
    expected_result = 1  # Only one locus meets all criteria

    # Run test
    result = calc_u(ref_gts, tgt_gts, [src_gts1, src_gts2], w, x, y_list)
    assert result == expected_result, f"Expected {expected_result}, got {result}"


def test_calc_q_with_two_sources():
    # Test data
    ref_gts = np.array([[0, 0, 0], [0, 1, 1], [1, 1, 1], [0, 0, 1]])
    tgt_gts = np.array([[0, 0, 0], [1, 1, 1], [1, 1, 1], [1, 1, 1]])
    src_gts1 = np.array([[0, 0, 0], [1, 1, 1], [1, 1, 1], [0, 0, 1]])
    src_gts2 = np.array([[1, 1, 1], [1, 1, 1], [0, 0, 0], [1, 1, 1]])
    w, y_list, quantile = 0.5, [1.0, 1.0], 0.95

    # Expected result: 95% quantile of the filtered tgt_gts frequencies
    expected_result = np.nan

    # Run test
    result = calc_q(ref_gts, tgt_gts, [src_gts1, src_gts2], w, y_list, quantile)
    if np.isnan(expected_result):
        assert np.isnan(result), f"Expected NaN, got {result}"
    else:
        assert np.isclose(
            result, expected_result
        ), f"Expected {expected_result}, got {result}"


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


def test_calc_seq_div():
    # Test case 1: Simple case with known divergence
    gts1 = np.array([[0, 1], [1, 0]])
    gts2 = np.array([[1, 0], [0, 1]])
    expected_divergence = np.array([[2, 0], [0, 2]])
    result = calc_seq_div(gts1, gts2)
    assert np.array_equal(
        result, expected_divergence
    ), f"Failed on test case 1 with result {result}"

    # Test case 2: Same populations (should result in zero divergence)
    gts1 = np.array([[1, 1], [1, 1]])
    gts2 = np.array([[1, 1], [1, 1]])
    expected_divergence = np.array([[0, 0], [0, 0]])
    result = calc_seq_div(gts1, gts2)
    assert np.array_equal(
        result, expected_divergence
    ), f"Failed on test case 2 with result {result}"

    # Test case 3:
    gts1 = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
    gts2 = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
    expected_divergence = np.array(
        [
            [0, 2, 2],
            [2, 0, 2],
            [2, 2, 0],
        ]
    )
    result = calc_seq_div(gts1, gts2)
    assert np.array_equal(
        result, expected_divergence
    ), f"Failed on test case 3 with result {result}"

    # Test case 4:
    gts1 = np.array(
        [
            [0, 1, 2],
            [1, 2, 0],
            [0, 2, 1],
        ]
    )
    gts2 = np.array(
        [
            [0, 1, 2],
            [1, 2, 0],
            [0, 2, 1],
        ]
    )
    expected_divergence = np.array(
        [
            [0, 3, 3],
            [3, 0, 3],
            [3, 3, 0],
        ]
    )
    result = calc_seq_div(gts1, gts2)
    assert np.array_equal(
        result, expected_divergence
    ), f"Failed on test case 4 with result {result}"


def test_calc_rd():
    # Test case 1
    src_gts = np.array([[0, 1], [1, 0]])
    ref_gts = np.array([[1, 0], [0, 1]])
    tgt_gts = np.array([[1, 1], [0, 0]])
    expected_ratio = 1 
    result = calc_rd(ref_gts, tgt_gts, src_gts)
    assert np.isclose(
        result, expected_ratio
    ), f"Failed on test case 1 with result {result}"

    # Test case 2
    src_gts = np.array([[0, 1], [1, 1]])
    ref_gts = np.array([[1, 0], [0, 1]])
    tgt_gts = np.array([[1, 0], [0, 1]])
    expected_ratio = 1.0
    result = calc_rd(ref_gts, tgt_gts, src_gts)
    assert np.isclose(
        result, expected_ratio
    ), f"Failed on test case 2 with result {result}"


def test_calc_abba_baba_basic():
    # Test data
    ref_gts = np.array([[0, 0, 1], [0, 0, 0], [1, 1, 1]])
    tgt_gts = np.array([[0, 1, 1], [0, 0, 1], [1, 1, 1]])
    src_gts = np.array([[1, 1, 1], [0, 1, 1], [1, 1, 1]])
    
    # Expected output
    expected_result = 1.0 

    # Run test
    result = calc_abba_baba(ref_gts, tgt_gts, src_gts)
    assert np.isclose(
        result, expected_result
    ), f"Expected {expected_result}, got {result}"


def test_calc_abba_baba_no_match():
    # Produces invalid output
    ref_gts = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    tgt_gts = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    src_gts = np.array([[1, 1, 1], [0, 1, 1], [1, 1, 1]])
    
    # Expected output
    expected_result = np.nan

    # Run test
    result = calc_abba_baba(ref_gts, tgt_gts, src_gts)
    if np.isnan(expected_result):
        assert np.isnan(result), f"Expected NaN, got {result}"
    else:
        assert np.isclose(
            result, expected_result
        ), f"Expected {expected_result}, got {result}"
