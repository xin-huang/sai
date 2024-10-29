# Copyright 2024 Xin Huang
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
from typing import Any
from sai.stats.features import calc_rd, calc_u, calc_q
from sai.utils.preprocessors import DataPreprocessor


class FeatureVectorsPreprocessor(DataPreprocessor):
    """
    A preprocessor subclass for generating feature vectors from genomic data.

    This class extends DataPreprocessor to include additional functionality for creating
    feature vectors based on genomic variants, reference and target individual genotypes,
    and window-based genomic statistics.
    """

    def __init__(
        self,
        w: float,
        x: float,
        y: list[float],
        output_file: str,
        quantile: float = 0.95,
    ):
        """
        Initializes FeatureVectorsPreprocessor with specific frequency thresholds
        and output file for storing generated feature vectors.

        Parameters
        ----------
        w : float
            Frequency threshold for `calc_u`.
        x : float
            Frequency threshold for `calc_u`.
        y : list[float]
            List of frequency thresholds for `calc_u` and `calc_q`.
        output_file : str
            Path to the output file to save processed feature vectors.
        """
        self.w = w
        self.x = x
        self.y = y
        self.quantile = quantile
        self.output_file = output_file

    def run(
        self,
        chr_name: str,
        ref_pop: str,
        tgt_pop: str,
        src_pop_list: list[str],
        start: int,
        end: int,
        ref_gts: np.ndarray,
        tgt_gts: np.ndarray,
        src_gts_list: list[np.ndarray],
        ploidy: int,
    ) -> dict[str, Any]:
        """
        Generates feature vectors for a specified genomic window.

        Parameters
        ----------
        chr_name : str
            Chromosome name.
        ref_pop : str
            Reference population name.
        tgt_pop : str
            Target population name.
        src_pop_list : list[str]
            List of source population names.
        start : int
            Start position of the genomic window.
        end : int
            End position of the genomic window.
        ref_gts : np.ndarray
            Genotype data for the reference population.
        tgt_gts : np.ndarray
            Genotype data for the target population.
        src_gts_list : list[np.ndarray]
            List of genotype arrays for each source population.
        ploidy: int
            Ploidy of the genome.

        Returns
        -------
        dict[str, Any]
            A dictionary containing calculated feature vectors for the genomic window.
        """
        items = {
            "chr_name": chr_name,
            "start": start,
            "end": end,
            "ref_pop": ref_pop,
            "tgt_pop": tgt_pop,
            "src_pop_list": src_pop_list,
        }

        # Calculate u statistic based on the provided w, x, and y
        items["u_statistic"] = calc_u(
            ref_gts=ref_gts,
            tgt_gts=tgt_gts,
            src_gts_list=src_gts_list,
            w=self.w,
            x=self.x,
            y_list=self.y,
            ploidy=ploidy,
        )

        # Calculate q statistic based on the provided w and y for each source population
        items["q_statistic"] = calc_q(
            ref_gts=ref_gts,
            tgt_gts=tgt_gts,
            src_gts_list=src_gts_list,
            w=self.w,
            y_list=self.y,
            quantile=self.quantile,
            ploidy=ploidy,
        )

        return items

    def process_items(self, items: dict[str, Any]) -> None:
        """
        Processes and writes a single dictionary of feature vectors to the output file.

        Parameters
        ----------
        items : dict[str, Any]
            A dictionary containing feature vectors for a genomic window.
        """
        with open(
            self.output_file, "a"
        ) as f:  # Open in append mode for continuous writing
            src_pop_str = ",".join(items["src_pop_list"])
            line = (
                f"{items['chr_name']}\t{items['start']}\t{items['end']}\t"
                f"{items['ref_pop']}\t{items['tgt_pop']}\t{src_pop_str}\t"
                f"{items['u_statistic']}\t{items['q_statistic']}\n"
            )
            f.write(line)
