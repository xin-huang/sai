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
from abc import ABC, abstractmethod
from typing import Dict, Any, Optional


class GenericStatistic(ABC):
    """
    Generic class for all statistics.

    This class provides a generic interface for implementing specific statistical measures
    from genotype matrices, typically representing different populations or samples.
    """

    def __init__(
        self,
        ref_gts: np.ndarray,
        tgt_gts: np.ndarray,
        ref_ploidy: int,
        tgt_ploidy: int,
        src_gts_list: list[np.ndarray],
        src_ploidy_list: list[int],
        out_gts: Optional[np.ndarray] = None,
        out_ploidy: Optional[int] = None,
    ):
        """
        Initializes the statistic with reference and target genotypes and their ploidies.

        Parameters
        ----------
        ref_gts : np.ndarray
            A 2D numpy array where each row represents a locus and each column represents an individual in the reference group.
        tgt_gts : np.ndarray
            A 2D numpy array where each row represents a locus and each column represents an individual in the target group.
        ref_ploidy : int
            Ploidy level of the reference population.
        tgt_ploidy : int
            Ploidy level of the target population.
        src_gts_list: list[np.ndarray]
            A list of 2D numpy arrays for each source population, where each row represents a locus and each column
            represents an individual in that source population.
        src_ploidy_list: list[int]
            A list of ploidy levels for the source populations. If provided, must match the number of source genotype arrays.
        out_gts: Optional[np.ndarray]
            A 2D numpy array where each row represents a locus and each column represents an individual in the outgroup.
            Default: None.
        out_ploidy: Optional[int]
            Ploidy level of the outgroup. Default: None.
        """
        self.ref_gts = ref_gts
        self.tgt_gts = tgt_gts
        self.src_gts_list = src_gts_list
        self.out_gts = out_gts
        self.ref_ploidy = ref_ploidy
        self.tgt_ploidy = tgt_ploidy
        self.src_ploidy_list = src_ploidy_list
        self.out_ploidy = out_ploidy

    @abstractmethod
    def compute(self, **kwargs) -> Dict[str, Any]:
        """
        Computes the statistic based on the input genotype data.

        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments specific to the statistic being implemented.

        Returns
        -------
        dict
            A dictionary containing the results of the statistic computation.
        """
        pass
