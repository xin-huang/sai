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
from sai.stats import DfStatistic


def test_DfStatistic_compute():
    ref_gts = np.array([[0, 0], [0, 0], [1, 1]])  # Referece population
    tgt_gts = np.array([[1, 0], [0, 1], [0, 1]])  # Taget population
    src_gts = np.array([[0, 1], [1, 0], [1, 0]])  # Source population
    # out_gts = None  # No outgroup provided

    # ref_freq = [0, 0, 1]
    # tgt_freq = [0.5, 0.5, 0.5]
    # src_freq = [0.5, 0.5, 0.5]
    # out_freq = [0, 0, 0]

    # pattern: 'abba'
    # site 0: (1-0)*0.5*0.5*(1-0) = 0.25
    # site 1: (1-0)*0.5*0.5*(1-0) = 0.25
    # site 2: (1-1)*0.5*0.5*(1-0) = 0
    # sum = 0.5

    # pattern: 'baba'
    # site 0: 0*(1-0.5)*0.5*(1-0) = 0
    # site 1: 0*(1-0.5)*0.5*(1-0) = 0
    # site 2: 1*(1-0.5)*0.5*(1-0) = 0.25
    # sum = 0.25

    # pattern: 'bbaa'
    # site 0: 0*0.5*(1-0.5)*(1-0) = 0
    # site 1: 0*0.5*(1-0.5)*(1-0) = 0
    # site 2: 1*0.5*(1-0.5)*(1-0) = 0.25
    # sum = 0.25

    # abba - baba = 0.5 - 0.25 = 0.25
    # abba + baba + 2 * bbaa = 0.5 + 0.25 + 2*0.25 = 1.25

    # Call the function with the test input
    df_stat = DfStatistic(
        ref_gts=ref_gts,
        tgt_gts=tgt_gts,
        src_gts_list=[src_gts],
        ref_ploidy=1,
        tgt_ploidy=1,
        src_ploidy_list=[1],
    )
    results = df_stat.compute()

    # Check the result
    expected_result = 0.2
    assert results["name"] == "df"
    assert np.isclose(
        results["value"][0], expected_result
    ), f"Expected {expected_result}, but got {results['value'][0]}"
