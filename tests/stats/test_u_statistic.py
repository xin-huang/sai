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
from sai.stats import UStatistic


def test_UStatistic_compute_basic():
    # Test data
    ref_gts = np.array([[0, 0, 1], [0, 0, 0], [1, 1, 1]])
    tgt_gts = np.array([[1, 1, 1], [1, 0, 0], [0, 1, 0]])
    src_gts = np.array([[0, 0, 0], [1, 1, 1], [1, 0, 1]])
    pos = np.array([0, 1, 2])
    w, x, y = 0.5, 0.5, ("=", 0)

    expected_result = 1  # Only the first site meets the criteria
    expected_positions = np.array([0])

    u_stat = UStatistic(
        ref_gts=ref_gts,
        tgt_gts=tgt_gts,
        src_gts_list=[src_gts],
        ref_ploidy=1,
        tgt_ploidy=1,
        src_ploidy_list=[1],
    )
    results = u_stat.compute(
        pos=pos,
        w=w,
        x=x,
        y_list=[y],
        anc_allele_available=False,
    )

    assert results["name"] == "U"
    assert (
        results["value"] == expected_result
    ), f"Expected {expected_result}, got {results['value']}"
    assert np.array_equal(results["cdd_pos"], expected_positions)

    w, x, y = 0.5, 0.5, ("=", 1)
    expected_result = 0
    expected_positions = np.array([])
    results = u_stat.compute(
        pos=pos,
        w=w,
        x=x,
        y_list=[y],
        anc_allele_available=True,
    )

    assert (
        results["value"] == expected_result
    ), f"Expected {expected_result}, got {results['value']}"
    assert np.array_equal(results["cdd_pos"], expected_positions)


def test_UStatistic_compute_no_match():
    # Test data with no matching sites
    ref_gts = np.array([[0, 1, 1], [1, 1, 1]])
    tgt_gts = np.array([[0, 0, 0], [1, 0, 1]])
    src_gts = np.array([[1, 1, 1], [1, 1, 1]])
    pos = np.array([0, 1])

    # Parameters for testing with no matching sites
    w, x, y = 0.3, 0.5, ("=", 0)

    expected_result = 0  # No sites meet the criteria
    expected_positions = np.array([])

    u_stat = UStatistic(
        ref_gts=ref_gts,
        tgt_gts=tgt_gts,
        src_gts_list=[src_gts],
        ref_ploidy=1,
        tgt_ploidy=1,
        src_ploidy_list=[1],
    )
    results = u_stat.compute(
        pos=pos,
        w=w,
        x=x,
        y_list=[y],
        anc_allele_available=False,
    )

    assert (
        results["value"] == expected_result
    ), f"Expected {expected_result}, got {results['value']}"
    assert np.array_equal(results["cdd_pos"], expected_positions)


def test_UStatistic_compute_all_match():
    # Test data where all sites meet the criteria
    ref_gts = np.array([[0, 0, 0], [0, 0, 0]])
    tgt_gts = np.array([[1, 1, 1], [1, 1, 1]])
    src_gts = np.array([[0, 0, 0], [0, 0, 0]])
    pos = np.array([0, 1])

    # Parameters for testing all sites matching
    w, x, y = 0.5, 0.5, ("=", 0)

    expected_result = 2  # All sites meet the criteria
    expected_positions = np.array([0, 1])

    u_stat = UStatistic(
        ref_gts=ref_gts,
        tgt_gts=tgt_gts,
        src_gts_list=[src_gts],
        ref_ploidy=1,
        tgt_ploidy=1,
        src_ploidy_list=[1],
    )
    results = u_stat.compute(
        pos=pos,
        w=w,
        x=x,
        y_list=[y],
        anc_allele_available=False,
    )

    assert (
        results["value"] == expected_result
    ), f"Expected {expected_result}, got {results['value']}"
    assert np.array_equal(results["cdd_pos"], expected_positions)


def test_UStatistic_compute_with_two_sources():
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

    u_stat = UStatistic(
        ref_gts=ref_gts,
        tgt_gts=tgt_gts,
        src_gts_list=[src_gts1, src_gts2],
        ref_ploidy=1,
        tgt_ploidy=1,
        src_ploidy_list=[1, 1],
    )
    results = u_stat.compute(
        pos=pos,
        w=w,
        x=x,
        y_list=y_list,
        anc_allele_available=False,
    )

    assert (
        results["value"] == expected_result
    ), f"Expected {expected_result}, got {results['value']}"
    assert np.array_equal(results["cdd_pos"], expected_positions)


def test_UStatistic_compute_with_mixed_ploidy():
    ref_gts = np.array([[0, 1, 0], [0, 1, 0], [2, 1, 0]])
    tgt_gts = np.array([[1, 1, 0], [1, 1, 1], [1, 1, 1]])
    src_gts = np.array([[0, 0, 0], [1, 1, 1], [0, 0, 0]])
    pos = np.array([0, 1, 2])
    w, x, y = 0.5, 0.5, ("=", 0)
    expected_result = 2
    expected_positions = np.array([0, 2])

    u_stat = UStatistic(
        ref_gts=ref_gts,
        tgt_gts=tgt_gts,
        src_gts_list=[src_gts],
        ref_ploidy=3,
        tgt_ploidy=1,
        src_ploidy_list=[2],
    )
    results = u_stat.compute(
        pos=pos,
        w=w,
        x=x,
        y_list=[y],
        anc_allele_available=False,
    )

    assert (
        results["value"] == expected_result
    ), f"Expected {expected_result}, got {results['value']}"
    assert np.array_equal(results["cdd_pos"], expected_positions)


def test_UStatistic_compute_with_missing_keys():
    ref_gts = np.array([[0, 1, 0], [0, 1, 0], [2, 1, 0]])
    tgt_gts = np.array([[1, 1, 0], [1, 1, 1], [1, 1, 1]])
    src_gts = np.array([[0, 0, 0], [1, 1, 1], [0, 0, 0]])
    pos = np.array([0, 1, 2])
    w, x, y = 0.5, 0.5, ("=", 0)

    with pytest.raises(ValueError):
        u_stat = UStatistic(
            ref_gts=ref_gts,
            tgt_gts=tgt_gts,
            src_gts_list=[src_gts],
            ref_ploidy=3,
            tgt_ploidy=1,
            src_ploidy_list=[2],
        )
        u_stat.compute(
            pos=pos,
            w=w,
            x=x,
            y_list=[y],
        )
