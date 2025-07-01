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
from sai.stats import FdStatistic


def test_FdStatistic_compute():
    ref_gts = np.array([[0, 1], [1, 0], [0, 1]])  # Reference population
    tgt_gts = np.array([[1, 0], [0, 1], [1, 0]])  # Target population
    src_gts = np.array([[1, 1], [1, 1], [1, 1]])  # Source population
    # out_gts = None  # No outgroup provided

    # ref_freq = [0.5, 0.5, 0.5]
    # tgt_freq = [0.5, 0.5, 0.5]
    # src_freq = [1, 1, 1]
    # out_freq = [0, 0, 0]

    # dnr_freq = src_freq = [1, 1, 1]
    # pattern: 'abba'
    # site 0: (1-0.5)*1*1*(1-0) = 0.5
    # site 1: (1-0.5)*1*1*(1-0) = 0.5
    # site 2: (1-0.5)*1*1*(1-0) = 0.5
    # sum = 1.5

    # pattern: 'baba'
    # site 0: 0.5*(1-1)*1*(1-0) = 0
    # site 1: 0.5*(1-1)*1*(1-0) = 0
    # site 2: 0.5*(1-1)*1*(1-0) = 0
    # sum = 0

    # abba_d - baba_d = 1.5

    # Call the function with the test input
    fd_stat = FdStatistic(
        ref_gts=ref_gts,
        tgt_gts=tgt_gts,
        src_gts_list=[src_gts],
        ref_ploidy=1,
        tgt_ploidy=1,
        src_ploidy_list=[1],
    )

    results = fd_stat.compute()

    # Check the result
    expected_result = 0
    assert results["name"] == "fd"
    assert np.isclose(
        results["value"][0], expected_result
    ), f"Expected {expected_result}, but got {results['value'][0]}"
