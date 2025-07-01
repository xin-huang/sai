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
from sai.stats.features import calc_u
from sai.stats.features import calc_q
from sai.stats.features import calc_freq
from sai.stats.features import compute_matching_loci


def test_calc_u_basic():
    # Test data
    ref_gts = np.array([[0, 0, 1], [0, 0, 0], [1, 1, 1]])
    tgt_gts = np.array([[1, 1, 1], [1, 0, 0], [0, 1, 0]])
    src_gts = np.array([[0, 0, 0], [1, 1, 1], [1, 0, 1]])
    pos = np.array([0, 1, 2])

    w, x, y = 0.5, 0.5, ("=", 0)
    expected_count = 1  # Only the first site meets the criteria
    expected_positions = np.array([0])
    result, loci_positions = calc_u(
        ref_gts, tgt_gts, [src_gts], pos, w, x, [y], [1, 1, 1]
    )
    assert result == expected_count
    assert np.array_equal(loci_positions, expected_positions)

    w, x, y = 0.5, 0.5, ("=", 1)
    expected_count = 0
    expected_positions = np.array([])
    result, loci_positions = calc_u(
        ref_gts, tgt_gts, [src_gts], pos, w, x, [y], [1, 1, 1], True
    )
    assert result == expected_count
    assert np.array_equal(loci_positions, expected_positions)


def test_calc_u_no_match():
    # Test data with no matching sites
    ref_gts = np.array([[0, 1, 1], [1, 1, 1]])
    tgt_gts = np.array([[0, 0, 0], [1, 0, 1]])
    src_gts = np.array([[1, 1, 1], [1, 1, 1]])
    pos = np.array([0, 1])

    # Parameters for testing with no matching sites
    w, x, y = 0.3, 0.5, ("=", 0)

    expected_count = 0  # No sites meet the criteria
    expected_positions = np.array([])
    result, loci_positions = calc_u(
        ref_gts, tgt_gts, [src_gts], pos, w, x, [y], [1, 1, 1]
    )
    assert result == expected_count
    assert np.array_equal(loci_positions, expected_positions)


def test_calc_u_all_match():
    # Test data where all sites meet the criteria
    ref_gts = np.array([[0, 0, 0], [0, 0, 0]])
    tgt_gts = np.array([[1, 1, 1], [1, 1, 1]])
    src_gts = np.array([[0, 0, 0], [0, 0, 0]])
    pos = np.array([0, 1])

    # Parameters for testing all sites matching
    w, x, y = 0.5, 0.5, ("=", 0)

    expected_count = 2  # All sites meet the criteria
    expected_positions = np.array([0, 1])

    result, loci_positions = calc_u(
        ref_gts, tgt_gts, [src_gts], pos, w, x, [y], [1, 1, 1]
    )
    assert result == expected_count
    assert np.array_equal(loci_positions, expected_positions)


def test_calc_q_basic():
    # Test data
    ref_gts = np.array([[0, 0, 1], [0, 0, 0], [1, 1, 1]])
    tgt_gts = np.array([[0, 1, 1], [0, 0, 1], [1, 1, 1]])
    src_gts = np.array([[1, 1, 1], [0, 1, 1], [1, 1, 1]])
    pos = np.array([0, 1, 2])
    w, y, quantile = 0.5, ("=", 1.0), 0.95

    # Expected output
    expected_result = 0.66667  # Only the first site meets the criteria
    expected_positions = np.array([0])

    # Run test
    result, loci_positions = calc_q(
        ref_gts, tgt_gts, [src_gts], pos, w, [y], [1, 1, 1], quantile
    )
    assert np.isclose(
        result, expected_result
    ), f"Expected {expected_result}, got {result}"
    assert np.array_equal(loci_positions, expected_positions)

    result, loci_positions = calc_q(
        ref_gts, tgt_gts, [src_gts], pos, w, [y], [1, 1, 1], quantile, True
    )
    assert np.isclose(
        result, expected_result
    ), f"Expected {expected_result}, got {result}"
    assert np.array_equal(loci_positions, expected_positions)


def test_calc_q_no_match():
    # Test data with no matching loci
    ref_gts = np.array([[0, 0, 1], [0, 0, 0]])
    tgt_gts = np.array([[0, 1, 1], [1, 1, 1]])
    src_gts = np.array([[1, 1, 1], [1, 1, 1]])
    pos = np.array([0, 1])
    w, y, quantile = (
        0.3,
        ("=", 0.0),
        0.95,
    )  # No tgt_gts frequencies < w and no src_gts frequencies == y

    # Expected output
    expected_positions = np.array([])

    # Run test
    result, loci_positions = calc_q(
        ref_gts, tgt_gts, [src_gts], pos, w, [y], [1, 1, 1], quantile
    )
    assert np.isnan(result), f"Expected NaN, got {result}"
    assert np.array_equal(loci_positions, expected_positions)


def test_calc_q_different_quantile():
    # Test data
    ref_gts = np.array([[0, 0, 1], [1, 0, 0], [0, 0, 1]])
    tgt_gts = np.array([[0, 1, 1], [1, 1, 1], [1, 1, 1]])
    src_gts = np.array([[0, 0, 0], [1, 1, 1], [1, 1, 1]])
    pos = np.array([0, 1, 2])
    w, y, quantile = 0.5, ("=", 1.0), 0.5

    # Expected output
    expected_result = (
        1.0  # 50% quantile (median) of [1.0, 1.0, 1.0] in tgt_gts that meets conditions
    )
    expected_positions = np.array([1, 2])

    # Run test
    result, loci_positions = calc_q(
        ref_gts, tgt_gts, [src_gts], pos, w, [y], [1, 1, 1], quantile
    )
    assert np.isclose(
        result, expected_result
    ), f"Expected {expected_result}, got {result}"
    assert np.array_equal(loci_positions, expected_positions)


def test_calc_q_edge_case():
    # Edge case where only one site meets criteria
    ref_gts = np.array([[0, 0, 1], [0, 0, 0], [1, 1, 1]])
    tgt_gts = np.array([[0, 1, 1], [1, 1, 1], [0, 0, 0]])
    src_gts = np.array([[0, 0, 0], [1, 1, 1], [1, 1, 1]])
    pos = np.array([0, 1, 2])
    w, y, quantile = 0.95, ("=", 1.0), 0.95

    # Expected output
    expected_result = 0.9666666666666667
    expected_positions = np.array([1])

    # Run test
    result, loci_positions = calc_q(
        ref_gts, tgt_gts, [src_gts], pos, w, [y], [1, 1, 1], quantile
    )
    assert np.isclose(
        result, expected_result
    ), f"Expected {expected_result}, got {result}"
    assert np.array_equal(loci_positions, expected_positions)


def test_calc_q_invalid_quantile():
    # Test data for invalid quantile
    ref_gts = np.array([[0, 0, 1]])
    tgt_gts = np.array([[0, 1, 1]])
    src_gts = np.array([[1, 1, 1]])
    pos = np.array([0])
    w, y, quantile = 0.5, ("=", 1.0), 1.5  # Invalid quantile (out of [0, 1] range)

    with pytest.raises(
        ValueError,
        match="Parameter quantile must be within the range \\[0, 1\\]",
    ):
        calc_q(ref_gts, tgt_gts, [src_gts], pos, w, [y], [1, 1, 1], quantile)


def test_calc_u_with_two_sources():
    # Test data
    ref_gts = np.array([[0, 0, 1], [0, 0, 0], [1, 1, 1]])
    tgt_gts = np.array([[0, 1, 1], [0, 0, 1], [1, 1, 1]])
    src_gts1 = np.array([[1, 1, 1], [0, 1, 1], [1, 1, 1]])
    src_gts2 = np.array([[1, 1, 1], [0, 1, 1], [1, 1, 1]])
    pos = np.array([0, 1, 2])
    w, x, y_list = 0.5, 0.5, [("=", 1.0), ("=", 1.0)]

    # Expected result: only loci where ref_freq < w, tgt_freq > x, src_freq1 == y1, and src_freq2 == y2
    expected_result = 1  # Only one locus meets all criteria
    expected_positions = np.array([0])

    # Run test
    result, loci_positions = calc_u(
        ref_gts, tgt_gts, [src_gts1, src_gts2], pos, w, x, y_list, [1, 1, 1]
    )
    assert result == expected_result, f"Expected {expected_result}, got {result}"
    assert np.array_equal(loci_positions, expected_positions)


def test_calc_q_with_two_sources():
    # Test data
    ref_gts = np.array([[1, 1, 0], [0, 1, 1], [1, 1, 1], [0, 0, 1]])
    tgt_gts = np.array([[0, 0, 0], [1, 1, 1], [1, 1, 1], [1, 1, 1]])
    src_gts1 = np.array([[0, 0, 0], [1, 1, 1], [1, 1, 1], [0, 0, 1]])
    src_gts2 = np.array([[1, 1, 1], [1, 1, 1], [0, 0, 0], [1, 1, 1]])
    pos = np.array([0, 1, 2, 3])
    w, y_list, quantile = 0.5, [("=", 1), ("=", 1)], 0.95

    # Expected result: 95% quantile of the filtered tgt_gts frequencies
    expected_result = np.nan
    expected_positions = np.array([])

    # Run test
    result, loci_positions = calc_q(
        ref_gts, tgt_gts, [src_gts1, src_gts2], pos, w, y_list, [1, 1, 1, 1], quantile
    )
    assert np.isnan(result), f"Expected NaN, got {result}"
    assert np.array_equal(loci_positions, expected_positions)

    w, x, y_list = 0.5, 0.8, [("=", 1), ("=", 0)]
    expected_result = 1
    expected_positions = np.array([0])

    result, loci_positions = calc_u(
        ref_gts, tgt_gts, [src_gts1, src_gts2], pos, w, x, y_list, [1, 1, 1, 1]
    )
    assert np.isclose(result, expected_result)
    assert np.array_equal(loci_positions, expected_positions)


def test_calc_u_mixed_ploidy():
    ref_gts = np.array([[0, 1, 0], [0, 1, 0], [2, 1, 0]])
    tgt_gts = np.array([[1, 1, 0], [1, 1, 1], [1, 1, 1]])
    src_gts = np.array([[0, 0, 0], [1, 1, 1], [0, 0, 0]])
    pos = np.array([0, 1, 2])
    w, x, y = 0.5, 0.5, ("=", 0)
    expected_count = 2
    expected_positions = np.array([0, 2])
    result, loci_positions = calc_u(
        ref_gts, tgt_gts, [src_gts], pos, w, x, [y], [3, 1, 2]
    )
    assert result == expected_count
    assert np.array_equal(loci_positions, expected_positions)


def test_calc_q_mixed_ploidy():
    ref_gts = np.array([[1, 1, 0], [0, 1, 1], [1, 1, 1], [0, 0, 1]])
    tgt_gts = np.array([[0, 0, 0], [1, 1, 1], [1, 1, 1], [1, 1, 1]])
    src_gts1 = np.array([[0, 0, 0], [1, 1, 1], [1, 1, 1], [0, 0, 1]])
    src_gts2 = np.array([[1, 1, 1], [1, 1, 1], [0, 0, 0], [1, 1, 1]])
    pos = np.array([0, 1, 2, 3])
    w, y_list, quantile = 0.5, [("=", 1), ("=", 1)], 0.95

    expected_result = np.nan
    expected_positions = np.array([])

    result, loci_positions = calc_q(
        ref_gts, tgt_gts, [src_gts1, src_gts2], pos, w, y_list, [2, 2, 4, 4], quantile
    )
    assert np.isnan(result), f"Expected NaN, got {result}"
    assert np.array_equal(loci_positions, expected_positions)


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
