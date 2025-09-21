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
from sai.stats import QStatistic


def test_QStatistic_compute_basic():
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
    q_stat = QStatistic(
        ref_gts=ref_gts,
        tgt_gts=tgt_gts,
        src_gts_list=[src_gts],
        ref_ploidy=1,
        tgt_ploidy=1,
        src_ploidy_list=[1],
    )
    results = q_stat.compute(
        pos=pos,
        w=w,
        y_list=[y],
        quantile=quantile,
        anc_allele_available=False,
    )

    assert results["name"] == "Q"
    assert np.isclose(
        results["value"], expected_result
    ), f"Expected {expected_result}, got {results['value']}"
    assert np.array_equal(results["cdd_pos"], expected_positions)

    results = q_stat.compute(
        pos=pos,
        w=w,
        y_list=[y],
        quantile=quantile,
        anc_allele_available=True,
    )

    assert np.isclose(
        results["value"], expected_result
    ), f"Expected {expected_result}, got {results['value']}"
    assert np.array_equal(results["cdd_pos"], expected_positions)


def test_QStatistic_compute_no_match():
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
    q_stat = QStatistic(
        ref_gts=ref_gts,
        tgt_gts=tgt_gts,
        src_gts_list=[src_gts],
        ref_ploidy=1,
        tgt_ploidy=1,
        src_ploidy_list=[1],
    )
    results = q_stat.compute(
        pos=pos,
        w=w,
        y_list=[y],
        quantile=quantile,
        anc_allele_available=False,
    )

    assert np.isnan(results["value"]), f"Expected NaN, got {results['value']}"
    assert np.array_equal(results["cdd_pos"], expected_positions)


def test_QStatistic_compute_different_quantile():
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
    q_stat = QStatistic(
        ref_gts=ref_gts,
        tgt_gts=tgt_gts,
        src_gts_list=[src_gts],
        ref_ploidy=1,
        tgt_ploidy=1,
        src_ploidy_list=[1],
    )
    results = q_stat.compute(
        pos=pos,
        w=w,
        y_list=[y],
        quantile=quantile,
        anc_allele_available=False,
    )

    assert np.isclose(
        results["value"], expected_result
    ), f"Expected {expected_result}, got {results['value']}"
    assert np.array_equal(results["cdd_pos"], expected_positions)


def test_QStatistic_compute_edge_case():
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
    q_stat = QStatistic(
        ref_gts=ref_gts,
        tgt_gts=tgt_gts,
        src_gts_list=[src_gts],
        ref_ploidy=1,
        tgt_ploidy=1,
        src_ploidy_list=[1],
    )
    results = q_stat.compute(
        pos=pos,
        w=w,
        y_list=[y],
        quantile=quantile,
        anc_allele_available=False,
    )

    assert np.isclose(
        results["value"], expected_result
    ), f"Expected {expected_result}, got {results['value']}"
    assert np.array_equal(results["cdd_pos"], expected_positions)


def test_QStatistic_compute_with_two_sources():
    # Test data
    ref_gts = np.array([[1, 1, 0], [0, 1, 1], [1, 1, 1], [0, 0, 1]])
    tgt_gts = np.array([[0, 0, 0], [1, 1, 1], [1, 1, 1], [1, 1, 1]])
    src_gts1 = np.array([[0, 0, 0], [1, 1, 1], [1, 1, 1], [0, 0, 1]])
    src_gts2 = np.array([[1, 1, 1], [1, 1, 1], [0, 0, 0], [1, 1, 1]])
    pos = np.array([0, 1, 2, 3])
    w, y_list, quantile = 0.5, [("=", 1), ("=", 1)], 0.95

    # Expected result: 95% quantile of the filtered tgt_gts frequencies
    expected_positions = np.array([])

    # Run test
    q_stat = QStatistic(
        ref_gts=ref_gts,
        tgt_gts=tgt_gts,
        src_gts_list=[src_gts1, src_gts2],
        ref_ploidy=1,
        tgt_ploidy=1,
        src_ploidy_list=[1, 1],
    )
    results = q_stat.compute(
        pos=pos,
        w=w,
        y_list=y_list,
        quantile=quantile,
        anc_allele_available=False,
    )

    assert np.isnan(results["value"]), f"Expected NaN, got {results['value']}"
    assert np.array_equal(results["cdd_pos"], expected_positions)


def test_QStatistic_compute_with_mixed_ploidy():
    ref_gts = np.array([[1, 1, 0], [0, 1, 1], [1, 1, 1], [0, 0, 1]])
    tgt_gts = np.array([[0, 0, 0], [1, 1, 1], [1, 1, 1], [1, 1, 1]])
    src_gts1 = np.array([[0, 0, 0], [1, 1, 1], [1, 1, 1], [0, 0, 1]])
    src_gts2 = np.array([[1, 1, 1], [1, 1, 1], [0, 0, 0], [1, 1, 1]])
    pos = np.array([0, 1, 2, 3])
    w, y_list, quantile = 0.5, [("=", 1), ("=", 1)], 0.95

    expected_positions = np.array([])

    q_stat = QStatistic(
        ref_gts=ref_gts,
        tgt_gts=tgt_gts,
        src_gts_list=[src_gts1, src_gts2],
        ref_ploidy=2,
        tgt_ploidy=2,
        src_ploidy_list=[4, 4],
    )
    results = q_stat.compute(
        pos=pos,
        w=w,
        y_list=y_list,
        quantile=quantile,
        anc_allele_available=False,
    )

    assert np.isnan(results["value"]), f"Expected NaN, got {results['value']}"
    assert np.array_equal(results["cdd_pos"], expected_positions)


def test_QStatistic_compute_with_missing_keys():
    ref_gts = np.array([[1, 1, 0], [0, 1, 1], [1, 1, 1], [0, 0, 1]])
    tgt_gts = np.array([[0, 0, 0], [1, 1, 1], [1, 1, 1], [1, 1, 1]])
    src_gts1 = np.array([[0, 0, 0], [1, 1, 1], [1, 1, 1], [0, 0, 1]])
    src_gts2 = np.array([[1, 1, 1], [1, 1, 1], [0, 0, 0], [1, 1, 1]])
    pos = np.array([0, 1, 2, 3])
    w, quantile = 0.5, 0.95

    with pytest.raises(ValueError):
        q_stat = QStatistic(
            ref_gts=ref_gts,
            tgt_gts=tgt_gts,
            src_gts_list=[src_gts1, src_gts2],
            ref_ploidy=2,
            tgt_ploidy=2,
            src_ploidy_list=[4, 4],
        )
        q_stat.compute(
            pos=pos,
            w=w,
            quantile=quantile,
            anc_allele_available=False,
        )
