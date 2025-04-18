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


from itertools import combinations, product
from typing import Iterator, Any
from sai.utils import read_data, split_genome
from sai.utils.generators import DataGenerator


class WindowGenerator(DataGenerator):
    """
    Generates genomic data for each specified window from VCF and other related files,
    allowing the user to select the number of source populations.
    """

    def __init__(
        self,
        vcf_file: str,
        chr_name: str,
        ref_ind_file: str,
        tgt_ind_file: str,
        src_ind_file: str,
        win_len: int,
        win_step: int,
        start: int = None,
        end: int = None,
        anc_allele_file: str = None,
        num_src: int = 1,
    ):
        """
        Initializes a new instance of WindowGenerator.

        Parameters
        ----------
        vcf_file : str
            The path to the VCF file containing variant data.
        chr_name: str
            The chromosome name to read from the VCF file.
        ref_ind_file : str
            The path to the file containing identifiers for reference populations.
        tgt_ind_file : str
            The path to the file containing identifiers for target populations.
        src_ind_file : str
            The path to the file containing identifiers for source populations.
        win_len : int
            The length of each window in base pairs.
        win_step : int
            The step size between windows in base pairs.
        start: int, optional
            The starting position (1-based, inclusive) on the chromosome. Default: None.
        end: int, optional
            The ending position (1-based, inclusive) on the chromosome. Default: None.
        anc_allele_file: str, optional
            Path to the file containing ancestral allele information. Default: None.
        num_src : int, optional
            The number of source populations to include in each combination. Default: 1.

        Raises
        ------
        ValueError
            If `win_len` is less than or equal to 0, if `win_step` is negative.
        """
        if win_len <= 0:
            raise ValueError("`win_len` must be greater than 0.")
        if win_step < 0:
            raise ValueError("`win_step` must be non-negative.")
        if num_src < 1:
            raise ValueError("`num_src` must be at least 1.")

        self.win_len = win_len
        self.win_step = win_step
        self.num_src = num_src
        self.chr_name = chr_name

        # Load data
        (
            self.ref_data,
            self.ref_samples,
            self.tgt_data,
            self.tgt_samples,
            self.src_data,
            self.src_samples,
            self.ploidy,
        ) = read_data(
            vcf_file=vcf_file,
            chr_name=self.chr_name,
            start=start,
            end=end,
            ref_ind_file=ref_ind_file,
            tgt_ind_file=tgt_ind_file,
            src_ind_file=src_ind_file,
            anc_allele_file=anc_allele_file,
            is_phased=False,
            filter_ref=False,
            filter_tgt=False,
            filter_src=False,
        )

        self.src_combinations = list(
            combinations(self.src_samples.keys(), self.num_src)
        )
        self.tgt_windows = {
            tgt_pop: split_genome(
                pos=(
                    self.tgt_data[tgt_pop].POS
                    if (start is None) and (end is None)
                    else [start, end - win_len + win_step]
                ),
                window_size=self.win_len,
                step_size=self.win_step,
                start=start,
            )
            for tgt_pop in self.tgt_samples
        }
        self.total_windows = sum(
            len(windows) * len(self.ref_samples) * len(self.src_combinations)
            for windows in self.tgt_windows.values()
        )

    def _window_generator(self) -> Iterator[dict[str, Any]]:
        """
        Generator function that yields genomic data for each window for each
        population combination, including specified source population combinations.

        Yields
        ------
        dict
            A dictionary containing population names, start and end positions,
            ploidy and phase information, reference, target, and source genotypes,
            and positions for each window.
        """
        for ref_pop, tgt_pop, src_comb in product(
            self.ref_samples, self.tgt_samples, self.src_combinations
        ):
            tgt_pos = self.tgt_data[tgt_pop].POS
            for start, end in self.tgt_windows[tgt_pop]:
                ref_gts = self.ref_data[ref_pop].GT[
                    (self.ref_data[ref_pop].POS >= start)
                    & (self.ref_data[ref_pop].POS <= end)
                ]
                tgt_gts = self.tgt_data[tgt_pop].GT[
                    (self.tgt_data[tgt_pop].POS >= start)
                    & (self.tgt_data[tgt_pop].POS <= end)
                ]
                src_gts_list = [
                    self.src_data[src_pop].GT[
                        (self.src_data[src_pop].POS >= start)
                        & (self.src_data[src_pop].POS <= end)
                    ]
                    for src_pop in src_comb
                ]

                sub_pos = tgt_pos[(tgt_pos >= start) & (tgt_pos <= end)]

                yield {
                    "chr_name": self.chr_name,
                    "ref_pop": ref_pop,
                    "tgt_pop": tgt_pop,
                    "src_pop_list": src_comb,  # List of source populations in this combination
                    "start": start,
                    "end": end,
                    "pos": sub_pos,
                    "ref_gts": ref_gts,
                    "tgt_gts": tgt_gts,
                    "src_gts_list": src_gts_list,  # List of genotypes for each source population in src_comb
                    "ploidy": self.ploidy,
                }

    def _none_window_generator(self) -> Iterator[dict[str, Any]]:
        """
        Generates empty window data when reference, target, or source data is missing.

        Yields
        ------
        dict[str, Any]
            A dictionary containing the following keys:
            - "chr_name" (str): The chromosome name.
            - "ref_pop" (str): Reference population name.
            - "tgt_pop" (str): Target population name.
            - "src_pop_list" (list[str]): List of source populations in this combination.
            - "start" (int): Start position of the window.
            - "end" (int): End position of the window.
            - "pos" (list[int]): Empty list, since no positions are available.
            - "ref_gts" (None): Placeholder for missing reference genotypes.
            - "tgt_gts" (None): Placeholder for missing target genotypes.
            - "src_gts_list" (None): Placeholder for missing source genotypes.
            - "ploidy" (None): Placeholder for missing ploidy information.
        """
        for ref_pop, tgt_pop, src_comb in product(
            self.ref_samples, self.tgt_samples, self.src_combinations
        ):
            for start, end in self.tgt_windows[tgt_pop]:
                yield {
                    "chr_name": self.chr_name,
                    "ref_pop": ref_pop,
                    "tgt_pop": tgt_pop,
                    "src_pop_list": src_comb,
                    "start": start,
                    "end": end,
                    "pos": [],
                    "ref_gts": None,
                    "tgt_gts": None,
                    "src_gts_list": None,
                    "ploidy": None,
                }

    def get(self) -> Iterator[dict[str, Any]]:
        """
        Returns the generator for window data.

        Returns
        -------
        generator
            A generator yielding genomic data for each window.
        """
        if (
            (self.ref_data is None)
            or (self.tgt_data is None)
            or (self.src_data is None)
        ):
            return self._none_window_generator()
        else:
            return self._window_generator()

    def __len__(self) -> int:
        """
        Returns the precomputed total number of windows across all population combinations.

        Returns
        -------
        int
            Total number of windows.
        """
        return self.total_windows
