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


@STAT_REGISTRY.register("U")
class UStatistic(GenericStatistic):
    """
    Class for computing the number of uniquely shared sites between the target and source populations (Racimo et al. 2017. Mol Biol Evol),
    conditional on allele frequency patterns in the reference and source populations.
    """

    STAT_NAME = "U"

    def compute(self, **kwargs) -> Dict[str, Any]:
        """
        Computes the count of genetic loci that meet specified allele frequency conditions
        across reference, target, and multiple source genotypes, with adjustments based on src_freq consistency.

        Parameters
        ----------
        pos : np.ndarray
            A 1D numpy array where each element represents the genomic position.
        w : float
            Threshold for the allele frequency in `ref_gts`. Only loci with frequencies less than `w` are counted.
            Must be within the range [0, 1].
        x : float
            Threshold for the allele frequency in `tgt_gts`. Only loci with frequencies greater than `x` are counted.
            Must be within the range [0, 1].
        y_list : list[float]
            List of exact allele frequency thresholds for each source population in `src_gts_list`.
            Must be within the range [0, 1] and have the same length as `src_gts_list`.
        anc_allele_available : bool
            If True, checks only for matches with `y` (assuming `1` represents the derived allele).
            If False, checks both matches with `y` and `1 - y`, taking the major allele in the source as the reference.

        Returns
        -------
        dict
            A dictionary containing:
            - 'name' : str
                The name of the statistic ("U").
            - 'value' : int
                The count of loci that meet all specified frequency conditions.
            - 'ccd_pos' : np.ndarray
                A 1D numpy array containing the genomic positions of the loci that meet the conditions.
        """
        required_keys = ["pos", "w", "x", "y_list", "anc_allele_available"]
        if missing := [k for k in required_keys if k not in kwargs]:
            raise ValueError(f"Missing required argument(s): {', '.join(missing)}")

        pos = kwargs["pos"]
        w = kwargs["w"]
        x = kwargs["x"]
        y_list = kwargs["y_list"]
        anc_allele_available = kwargs["anc_allele_available"]
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

        # Apply final conditions
        condition &= tgt_freq > x

        loci_indices = np.where(condition)[0]
        loci_positions = pos[loci_indices]
        count = loci_indices.size

        # Return count of matching loci
        return {"name": self.STAT_NAME, "value": count, "cdd_pos": loci_positions}
