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


import numpy as np
from sai.stats import DancStatistic


def test_DancStatistic_compute():
    ref_gts = np.array([[0, 0], [0, 0], [1, 1]])  # Referece population
    tgt_gts = np.array([[1, 0], [0, 1], [0, 1]])  # Taget population
    src_gts = np.array([[0, 1], [1, 0], [1, 0]])  # Source population
    # out_gts = None  # No outgroup provided

    # D ancestral
    # baaa - abaa = 0.25 - 0.5 = -0.25
    # baaa + abaa = 0.25 + 0.5 = 0.75
    # (baaa - abaa) / (baaa + abaa) = -1/3

    danc_stat = DancStatistic(
        ref_gts=ref_gts,
        tgt_gts=tgt_gts,
        src_gts_list=[src_gts],
        ref_ploidy=1,
        tgt_ploidy=1,
        src_ploidy_list=[1],
    )
    results = danc_stat.compute()

    expected_result = -1 / 3

    assert results["name"] == "Danc"
    assert np.isclose(
        results["value"][0], expected_result
    ), f"Expected {expected_result}, but got {results['value'][0]}"
