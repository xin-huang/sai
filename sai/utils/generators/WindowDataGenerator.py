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


from itertools import combinations, product
from sai.utils import read_data, split_genome
from sai.utils.generators import DataGenerator


class WindowDataGenerator(DataGenerator):
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
        num_src: int = 1,
        anc_allele_file: str = None,
        ploidy: int = 2,
        is_phased: bool = True,
    ):
        """
        Initializes a new instance of WindowDataGenerator.

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
        num_src : int, optional
            The number of source populations to include in each combination. Default: 1.
        anc_allele_file : str, optional
            The path to the file containing ancestral allele information.
            Default: None.
        ploidy : int, optional
            The ploidy of the genome. Default: 2.
        is_phased : bool, optional
            Specifies whether the genotype data is phased. Default: True.

        Raises
        ------
        ValueError
            If `win_len` is less than or equal to 0, if `win_step` is negative,
            or if `ploidy` is less than or equal to 0.
        """
        if win_len <= 0:
            raise ValueError("`win_len` must be greater than 0.")
        if win_step < 0:
            raise ValueError("`win_step` must be non-negative.")
        if ploidy <= 0:
            raise ValueError("`ploidy` must be greater than 0.")
        if num_src < 1:
            raise ValueError("`num_src` must be at least 1.")

        if is_phased:
            self.ploidy = 1
        else:
            self.ploidy = ploidy
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
        ) = read_data(
            vcf_file=vcf_file,
            chr_name=self.chr_name,
            ref_ind_file=ref_ind_file,
            tgt_ind_file=tgt_ind_file,
            src_ind_file=src_ind_file,
            anc_allele_file=anc_allele_file,
            is_phased=is_phased,
            filter_ref=False,
            filter_tgt=False,
            filter_src=False,
        )

        # Precompute the total number of windows across all combinations
        self.total_windows = 0
        for _, tgt_pop, src_comb in product(
            self.ref_data,
            self.tgt_data,
            combinations(self.src_data.keys(), self.num_src),
        ):
            tgt_pos = self.tgt_data[tgt_pop].POS
            windows = split_genome(
                pos=tgt_pos,
                window_size=self.win_len,
                step_size=self.win_step,
            )
            self.total_windows += len(windows)

    def window_data_generator(self):
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
            self.ref_data,
            self.tgt_data,
            combinations(self.src_data.keys(), self.num_src),
        ):
            tgt_pos = self.tgt_data[tgt_pop].POS

            # Split the genome data into windows for each position
            windows = split_genome(
                pos=tgt_pos,
                window_size=self.win_len,
                step_size=self.win_step,
            )

            for start, end in windows:
                ref_gts = self.ref_data[ref_pop].GT[
                    (self.ref_data[ref_pop].POS > start)
                    & (self.ref_data[ref_pop].POS <= end)
                ]
                tgt_gts = self.tgt_data[tgt_pop].GT[
                    (tgt_pos > start) & (tgt_pos <= end)
                ]

                # Get genotypes for each source population in the combination
                src_gts_list = [
                    self.src_data[src_pop].GT[
                        (self.src_data[src_pop].POS > start)
                        & (self.src_data[src_pop].POS <= end)
                    ]
                    for src_pop in src_comb
                ]

                sub_pos = tgt_pos[(tgt_pos > start) & (tgt_pos <= end)]

                yield {
                    "chr_name": self.chr_name,
                    "ref_pop": ref_pop,
                    "tgt_pop": tgt_pop,
                    "src_pop_list": src_comb,  # List of source populations in this combination
                    "start": start,
                    "end": end,
                    "ref_gts": ref_gts,
                    "tgt_gts": tgt_gts,
                    "src_gts_list": src_gts_list,  # List of genotypes for each source population in src_comb
                    "ploidy": self.ploidy,
                }

    def get(self):
        """
        Returns the generator for window data.

        Returns
        -------
        generator
            A generator yielding genomic data for each window.
        """
        return self.window_data_generator()

    def __len__(self):
        """
        Returns the precomputed total number of windows across all population combinations.

        Returns
        -------
        int
            Total number of windows.
        """
        return self.total_windows
