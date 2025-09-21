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
from sai.stats.stat_utils import compute_matching_loci


@STAT_REGISTRY.register("Q")
class QStatistic(GenericStatistic):
    """
    Class for computing the quantile statistic in the target population (Racimo et al. 2017. Mol Biol Evol),
    conditional on allele frequency patterns in the reference and source populations.
    """

    STAT_NAME = "Q"

    def compute(self, **kwargs) -> Dict[str, Any]:
        """
        Calculates a specified quantile of derived allele frequencies in `tgt_gts` for loci that meet specific conditions
        across reference and multiple source genotypes, with adjustments based on src_freq consistency.

        Parameters
        ----------
        pos: np.ndarray
            A 1D numpy array where each element represents the genomic position.
        w : float
            Frequency threshold for the derived allele in `ref_gts`. Only loci with frequencies lower than `w` are included.
            Must be within the range [0, 1].
        y_list : list[float]
            List of exact frequency thresholds for each source population in `src_gts_list`.
            Must be within the range [0, 1] and have the same length as `src_gts_list`.
        quantile : float
            The quantile to compute for the filtered `tgt_gts` frequencies. Must be within the range [0, 1].
        anc_allele_available : bool
            If True, checks only for matches with `y` (assuming `1` represents the derived allele).
            If False, checks both matches with `y` and `1 - y`, taking the major allele in the source as the reference.

        Returns
        -------
        dict
            A dictionary containing:
            - 'name' : str
                The name of the statistic ("Q").
            - 'value' : float
                The specified quantile of the derived allele frequencies in `tgt_gts` for loci meeting the specified conditions,
                or NaN if no loci meet the criteria.
            - 'ccd_pos' : np.ndarray
                A 1D numpy array containing the genomic positions of the loci that meet the conditions.
        """
        required_keys = ["pos", "w", "y_list", "anc_allele_available", "quantile"]
        if missing := [k for k in required_keys if k not in kwargs]:
            raise ValueError(f"Missing required argument(s): {', '.join(missing)}")

        pos = kwargs["pos"]
        w = kwargs["w"]
        y_list = kwargs["y_list"]
        anc_allele_available = kwargs["anc_allele_available"]
        quantile = kwargs["quantile"]
        ploidy = [self.ref_ploidy, self.tgt_ploidy] + self.src_ploidy_list

        ref_freq, tgt_freq, condition = compute_matching_loci(
            self.ref_gts,
            self.tgt_gts,
            self.src_gts_list,
            w,
            y_list,
            ploidy,
            anc_allele_available,
        )

        # Filter `tgt_gts` frequencies based on the combined condition
        filtered_tgt_freq = tgt_freq[condition]
        filtered_positions = pos[condition]

        # Return NaN if no loci meet the criteria
        if filtered_tgt_freq.size == 0:
            threshold = np.nan
            loci_positions = np.array([])
        else:
            threshold = np.nanquantile(filtered_tgt_freq, quantile)
            loci_positions = filtered_positions[filtered_tgt_freq >= threshold]

        # Calculate and return the specified quantile of the filtered `tgt_gts` frequencies
        return {"name": self.STAT_NAME, "value": threshold, "cdd_pos": loci_positions}
