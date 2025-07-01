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
from scipy.spatial.distance import cdist
from typing import Dict, Any
from sai.registries.stat_registry import STAT_REGISTRY
from sai.stats import GenericStatistic


@STAT_REGISTRY.register("DD")
class DdStatistic(GenericStatistic):
    """
    Class for computing the average difference of the sequence divergence.

    The DD statistic quantifies the difference in average pairwise sequence
    divergence between a source population and two target populations (reference
    and target), using Manhattan (cityblock) distance.
    """

    STAT_NAME = "DD"

    def compute(self, **kwargs) -> Dict[str, Any]:
        """
        Computes the DD statistic for each source population.

        For each source population, the method calculates pairwise Manhattan distances
        between the source and both the target and reference populations, averages the
        distances per genome, and computes the difference in mean divergence.

        Parameters
        ----------
        **kwargs : dict
            Unused. Present to maintain compatibility with the base class interface.

        Returns
        -------
        dict
            A dictionary containing:
            - 'name' : str
                The name of the statistic ("DD").
            - 'value' : list[float]
                A list of DD values, one for each source population.
        """
        dd_results = []

        for i in range(len(self.src_gts_list)):
            # pairwise distances
            src_gts = self.src_gts_list[i]
            seq_divs_src_tgt = cdist(src_gts.T, self.tgt_gts.T, metric="cityblock")
            seq_divs_src_ref = cdist(src_gts.T, self.ref_gts.T, metric="cityblock")

            # mean of each row
            mean_src_tgt = np.mean(seq_divs_src_tgt, axis=1)
            mean_src_ref = np.mean(seq_divs_src_ref, axis=1)

            dd = np.mean(mean_src_ref - mean_src_tgt)
            dd_results.append(dd)

        return {"name": self.STAT_NAME, "value": dd_results}
