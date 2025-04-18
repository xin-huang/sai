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


import pysam
from typing import Iterator
from sai.utils import split_genome
from sai.utils.generators import DataGenerator


class ChunkGenerator(DataGenerator):
    """
    Generates genome chunks from VCF windows for parallel processing.

    This class splits genomic windows into non-overlapping chunks assigned to workers,
    based on the VCF file length and a user-defined window and step size.
    """

    def __init__(
        self,
        vcf_file: str,
        chr_name: str,
        step_size: int,
        window_size: int,
        num_chunks: int,
    ):
        """
        Initializes a new instance of ChunkGenerator.

        Parameters
        ----------
        vcf_file : str
            Path to the VCF file to process.
        chr_name: str
            Name of the chromosome to process.
        step_size : int
            Step size for generating windows.
        window_size : int
            Window size for generating windows.
        num_chunks : int
            Number of chunks to split the windows among.

        Raises
        ------
        ValueError
            If the specified chromosome is not found in the VCF file.
        """
        with pysam.VariantFile(vcf_file) as vcf:
            first_pos = last_pos = None
            for rec in vcf:
                if rec.chrom != chr_name:
                    if first_pos is not None:
                        break
                    continue
                if first_pos is None:
                    first_pos = rec.pos
                last_pos = rec.pos

        if first_pos is None:
            raise ValueError(f"Chromosome {chr_name} not found in VCF.")

        windows = split_genome([first_pos, last_pos], window_size, step_size)

        self.chunks = self._split_windows_ranges(windows, num_chunks)
        self.num_chunks = len(self.chunks)
        self.chr_name = chr_name

    def get(self) -> Iterator[tuple[str, int, int]]:
        """
        Yields a tuple representing the chunk assigned to each worker.

        Yields
        ------
        tuple of int
            A tuple representing the range (chr_name, start, end) assigned to each worker.
        """
        for chunk in self.chunks:
            yield {
                "chr_name": self.chr_name,
                "start": chunk[0],
                "end": chunk[1],
            }

    def __len__(self) -> int:
        """
        Returns the number of chunks.

        Returns
        -------
        int
            Number of chunks.
        """
        return self.num_chunks

    def _split_windows_ranges(self, windows: list, num_chunks: int) -> list:
        """
        Splits the list of windows into ranges assigned to each chunk.

        Each range is defined by the first window's start and the last window's end
        within that split.

        Parameters
        ----------
        windows : list of tuple
            List of (start, end) tuples representing windows.
        num_chunks : int
            Number of chunks to divide the windows among.

        Returns
        -------
        list of tuple
            List of (start, end) tuples representing the ranges for each chunk.
        """
        avg = len(windows) // num_chunks
        remainder = len(windows) % num_chunks
        result = []
        start_idx = 0

        for i in range(num_chunks):
            end_idx = start_idx + avg + (1 if i < remainder else 0)
            sub = windows[start_idx:end_idx]
            if sub:
                result.append((sub[0][0], sub[-1][1]))
            start_idx = end_idx

        return result
