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


from typing import Any
from sai.utils.generators import WindowGenerator
from sai.utils.preprocessors import DataPreprocessor
from .feature_preprocessor import FeaturePreprocessor


class ChunkPreprocessor(DataPreprocessor):
    """
    Preprocesses VCF data in genomic windows and applies feature preprocessing.

    This class generates genomic windows from a VCF file, processes them
    with specified reference, target, and source individuals, and computes
    feature vectors using the provided feature preprocessor.
    """

    def __init__(
        self,
        vcf_file: str,
        ref_ind_file: str,
        tgt_ind_file: str,
        src_ind_file: str,
        win_len: int,
        win_step: int,
        w: float,
        y: list[float],
        output_file: str,
        stat_type: str,
        anc_allele_file: str = None,
        num_src: int = 1,
    ):
        """
        Initializes a new instance of ChunkPreprocessor.

        Parameters
        ----------
        vcf_file : str
            Path to the VCF file to process.
        ref_ind_file : str
            Path to the file containing reference individual IDs.
        tgt_ind_file : str
            Path to the file containing target individual IDs.
        src_ind_file : str
            Path to the file containing source individual IDs.
        win_len : int
            Window length for generating genomic windows.
        win_step : int
            Step size for sliding windows across the genome.
        w : float
            Parameter w for feature vector computation.
        y : list of float
            List of y parameters for feature vector computation.
        output_file : str
            Path to the output file for storing feature vectors.
        stat_type : str
            Type of statistic to compute for feature vectors.
        anc_allele_file : str, optional
            Path to the ancestral allele file. If None, ancestral allele
            information is considered unavailable.
        num_src : int, optional
            Number of source populations to use. Default is 1.
        """
        self.vcf_file = vcf_file
        self.ref_ind_file = ref_ind_file
        self.tgt_ind_file = tgt_ind_file
        self.src_ind_file = src_ind_file
        self.win_len = win_len
        self.win_step = win_step
        self.anc_allele_file = anc_allele_file
        self.num_src = num_src

        anc_allele_available = anc_allele_file is not None

        self.feature_preprocessor = FeaturePreprocessor(
            w=w,
            y=y,
            output_file=output_file,
            stat_type=stat_type,
            anc_allele_available=anc_allele_available,
        )

    def run(self, chr_name: str, start: int, end: int) -> list[dict[str, Any]]:
        """
        Runs the preprocessing pipeline on a specific chromosome region.

        Generates genomic windows within the specified chromosome region,
        processes each window to compute feature vectors, and aggregates the results.

        Parameters
        ----------
        chr_name : str
            Name of the chromosome to process.
        start : int
            Start position (1-based, inclusive) of the region to process.
        end : int
            End position (1-based, exclusive) of the region to process.

        Returns
        -------
        list of dict of {str: Any}
            A list of dictionaries containing computed feature vectors for each genomic window.
        """
        window_generator = WindowGenerator(
            vcf_file=self.vcf_file,
            chr_name=chr_name,
            start=start,
            end=end,
            ref_ind_file=self.ref_ind_file,
            tgt_ind_file=self.tgt_ind_file,
            src_ind_file=self.src_ind_file,
            win_len=self.win_len,
            win_step=self.win_step,
            anc_allele_file=self.anc_allele_file,
            num_src=self.num_src,
        )

        items = []

        for item in window_generator.get():
            items.extend(self.feature_preprocessor.run(**item))

        return items

    def process_items(self, items: list[dict[str, Any]]) -> None:
        """
        Processes and writes computed feature vectors to the output.

        Parameters
        ----------
        items : list of dict of {str: Any}
            A list of dictionaries containing computed feature vectors for each genomic window.
        """
        self.feature_preprocessor.process_items(items)
