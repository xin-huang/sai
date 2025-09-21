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
from typing import Dict, Any
from sai.registries.stat_registry import STAT_REGISTRY
from sai.stats import GenericStatistic
from sai.stats.stat_utils import calc_four_pops_freq, calc_pattern_sum


@STAT_REGISTRY.register("df")
class DfStatistic(GenericStatistic):
    """
    Class for computing the distance fraction (df) statistic (Pfeifer and Kapan. 2019. BMC Bioinformatics).

    The df statistic quantifies the relative excess of shared derived alleles
    using ABBA, BABA, and BBAA site patterns across a four-population test.
    """

    STAT_NAME = "df"

    def compute(self, **kwargs) -> Dict[str, Any]:
        """
        Computes the df statistic for each source population.

        This method computes df for each source population using site pattern
        counts based on allele frequency input.

        Parameters
        ----------
        **kwargs : dict
            Unused. Present to maintain compatibility with the base class interface.

        Returns
        -------
        dict
            A dictionary containing:
            - 'name' : str
                The name of the statistic ("df").
            - 'value' : list[float]
                A list of df values, one per source population.
        """
        df_results = []

        for i in range(len(self.src_gts_list)):
            ref_freq, tgt_freq, src_freq, out_freq = calc_four_pops_freq(
                ref_gts=self.ref_gts,
                tgt_gts=self.tgt_gts,
                src_gts=self.src_gts_list[i],
                out_gts=self.out_gts,
                ref_ploidy=self.ref_ploidy,
                tgt_ploidy=self.tgt_ploidy,
                src_ploidy=self.src_ploidy_list[i],
                out_ploidy=self.out_ploidy,
            )

            abba = calc_pattern_sum(ref_freq, tgt_freq, src_freq, out_freq, "abba")
            baba = calc_pattern_sum(ref_freq, tgt_freq, src_freq, out_freq, "baba")
            bbaa = calc_pattern_sum(ref_freq, tgt_freq, src_freq, out_freq, "bbaa")

            numerator = abba - baba
            denominator = abba + baba + 2 * bbaa

            df = numerator / denominator if denominator != 0 else np.nan
            df_results.append(df)

        return {"name": self.STAT_NAME, "value": df_results}
