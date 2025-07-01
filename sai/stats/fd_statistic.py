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


@STAT_REGISTRY.register("fd")
class FdStatistic(GenericStatistic):
    """
    Class for computing the dynamic estimator of the proportion of introgression (Martin et al. 2015. Mol Biol Evol).

    The fd statistic is a dynamic estimator of the proportion of introgression
    from a source into the target population, based on ABBA-BABA pattern sums.
    """

    STAT_NAME = "fd"

    def compute(self, **kwargs) -> Dict[str, Any]:
        """
        Computes the fd statistic for each source population.

        This method iterates over each source population, computes allele
        frequencies, and evaluates the fd value using the standard numerator and
        dynamic denominator formulation.

        Parameters
        ----------
        **kwargs : dict
            Unused. Present to maintain compatibility with the base class interface.

        Returns
        -------
        dict
            A dictionary containing:
            - 'name' : str
                The name of the statistic ("fd").
            - 'value' : list[float]
                List of fd values, one for each source population.
        """
        fd_results = []

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

            abba_n = calc_pattern_sum(ref_freq, tgt_freq, src_freq, out_freq, "abba")
            baba_n = calc_pattern_sum(ref_freq, tgt_freq, src_freq, out_freq, "baba")

            dnr_freq = np.maximum(tgt_freq, src_freq)

            abba_d = calc_pattern_sum(ref_freq, dnr_freq, dnr_freq, out_freq, "abba")
            baba_d = calc_pattern_sum(ref_freq, dnr_freq, dnr_freq, out_freq, "baba")

            numerator = abba_n - baba_n
            denominator = abba_d - baba_d

            fd = numerator / denominator if denominator != 0 else np.nan
            fd_results.append(fd)

            # for i in range(len(ref_freq)):
            #    print(f"{ref_freq[i]}\t{tgt_freq[i]}\t{src_freq[i]}\t{out_freq[i]}")

        return {"name": self.STAT_NAME, "value": fd_results}
