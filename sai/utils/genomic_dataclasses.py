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


import allel
from dataclasses import dataclass


@dataclass
class ChromosomeData:
    """
    A data structure for storing chromosome-specific genotype information.

    Attributes
    ----------
    REF : list[str]
        A list of reference alleles for each variant position.
    ALT : list[str]
        A list of alternate alleles for each variant position.
    POS : list[int]
        A list of genomic positions corresponding to each variant.
    GT : list[allel.GenotypeVector]
        A list of genotype vectors, where each vector represents the genotype
        information for a specific variant position.
    """

    REF: list[str]
    ALT: list[str]
    POS: list[int]
    GT: list[allel.GenotypeVector]


@dataclass
class GenomicSegment:
    """
    Represents a segment of genomic data with relevant genotype and positional information.

    Parameters
    ----------
    chr_name : str
        Chromosome name of the genomic segment.
    start : int
        Start position of the genomic segment.
    end : int
        End position of the genomic segment.
    ploidy : int
        Ploidy level for this genomic segment (e.g., 2 for diploid).
    is_phased : bool
        Indicates if the segment is phased.
    ref_gts : list[allel.GenotypeVector]
        List of genotype vectors for the reference population across the segment.
    tgt_gts : list[allel.GenotypeVector]
        List of genotype vectors for the target population across the segment.
    src_gts : list[allel.GenotypeVector]
        List of genotype vectors from the source population across the segment.
    pos : List[int]
        List of genomic positions within the segment.
    """

    chr_name: str
    start: int
    end: int
    ploidy: int
    is_phased: bool
    ref_gts: list[allel.GenotypeVector]
    tgt_gts: list[allel.GenotypeVector]
    src_gts: list[allel.GenotypeVector]
    pos: list[int]
