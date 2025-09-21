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
from typing import Tuple, Optional, Union


def calc_freq(gts: np.ndarray, ploidy: int = 1) -> np.ndarray:
    """
    Calculates allele frequencies, supporting both phased and unphased data.

    Parameters
    ----------
    gts : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual.
    ploidy : int, optional
        Ploidy level of the organism. If ploidy=1, the function assumes phased data and calculates
        frequency by taking the mean across individuals. For unphased data, it calculates frequency by
        dividing the sum across individuals by the total number of alleles. Default is 1.

    Returns
    -------
    np.ndarray
        An array of allele frequencies for each locus.

    Raises
    ------
    ValueError
        If ploidy is not a positive integer.
    """
    if not isinstance(ploidy, int) or ploidy <= 0:
        raise ValueError("ploidy must be a positive integer.")

    return np.sum(gts, axis=1) / (gts.shape[1] * ploidy)


def compute_matching_loci(
    ref_gts: np.ndarray,
    tgt_gts: np.ndarray,
    src_gts_list: list[np.ndarray],
    w: float,
    y_list: list[tuple[str, float]],
    ploidy: list[int],
    anc_allele_available: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Computes loci that meet specified allele frequency conditions across reference, target, and source genotypes.

    Parameters
    ----------
    ref_gts : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual in the reference group.
    tgt_gts : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual in the target group.
    src_gts_list : list of np.ndarray
        A list of 2D numpy arrays for each source population, where each row represents a locus and each column
        represents an individual in that source population.
    w : float
        Threshold for the allele frequency in `ref_gts`. Only loci with frequencies less than `w` are counted.
        Must be within the range [0, 1].
    y_list : list of tuple[str, float]
        List of allele frequency conditions for each source population in `src_gts_list`.
        Each entry is a tuple (operator, threshold), where:
        - `operator` can be '=', '<', '>', '<=', '>='
        - `threshold` is a float within [0, 1]
        The length must match `src_gts_list`.
    ploidy : list[int]
        Ploidy values for reference, target, and one or more source populations (in that order).
    anc_allele_available : bool
        If True, checks only for matches with `y` (assuming `1` represents the derived allele).
        If False, checks both matches with `y` and `1 - y`, taking the dominant allele in the source as the reference.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        - Adjusted reference allele frequencies (`ref_freq`).
        - Adjusted target allele frequencies (`tgt_freq`).
        - Boolean array indicating loci that meet the specified frequency conditions (`condition`).
    """
    # Validate input parameters
    if not (0 <= w <= 1):
        raise ValueError("Parameters w must be within the range [0, 1].")

    for op, y in y_list:
        if not (0 <= y <= 1):
            raise ValueError(f"Invalid value in y_list: {y}. within the range [0, 1].")
        if op not in ("=", "<", ">", "<=", ">="):
            raise ValueError(
                f"Invalid operator in y_list: {op}. Must be '=', '<', '>', '<=', or '>='."
            )

    if len(src_gts_list) != len(y_list):
        raise ValueError("The length of src_gts_list and y_list must match.")

    # Compute allele frequencies
    ref_freq = calc_freq(ref_gts, ploidy[0])
    tgt_freq = calc_freq(tgt_gts, ploidy[1])
    src_freq_list = [
        calc_freq(src_gts, ploidy_val)
        for src_gts, ploidy_val in zip(src_gts_list, ploidy[2:])
    ]

    # Check match for each `y`
    op_funcs = {
        "=": lambda src_freq, y: src_freq == y,
        "<": lambda src_freq, y: src_freq < y,
        ">": lambda src_freq, y: src_freq > y,
        "<=": lambda src_freq, y: src_freq <= y,
        ">=": lambda src_freq, y: src_freq >= y,
    }

    match_conditions = [
        op_funcs[op](src_freq, y) for src_freq, (op, y) in zip(src_freq_list, y_list)
    ]
    all_match_y = np.all(match_conditions, axis=0)

    if not anc_allele_available:
        # Check if all source populations match `1 - y`
        match_conditions_1_minus_y = [
            op_funcs[op](src_freq, 1 - y)
            for src_freq, (op, y) in zip(src_freq_list, y_list)
        ]
        all_match_1_minus_y = np.all(match_conditions_1_minus_y, axis=0)
        all_match = all_match_y | all_match_1_minus_y

        # Identify loci where all sources match `1 - y` for frequency inversion
        inverted = all_match_1_minus_y

        # Invert frequencies for these loci
        ref_freq[inverted] = 1 - ref_freq[inverted]
        tgt_freq[inverted] = 1 - tgt_freq[inverted]
    else:
        all_match = all_match_y

    # Final condition: locus must satisfy source matching and have `ref_freq < w`
    condition = all_match & (ref_freq < w)

    return ref_freq, tgt_freq, condition


def calc_four_pops_freq(
    ref_gts: np.ndarray,
    tgt_gts: np.ndarray,
    src_gts: np.ndarray,
    out_gts: Optional[np.ndarray] = None,
    ref_ploidy: int = 1,
    tgt_ploidy: int = 1,
    src_ploidy: int = 1,
    out_ploidy: int = 1,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculates allele frequencies for four populations given their genotype matrices.

    Parameters
    ----------
    ref_gts : np.ndarray
        Genotype matrix for the reference population.
    tgt_gts : np.ndarray
        Genotype matrix for the target population.
    src_gts : np.ndarray
        Genotype matrix for the source population.
    out_gts : np.ndarray
        Genotype matrix for the outgroup. If None, the outgroup frequency is assumed to be 0 at all loci.
        Default: None.
    ref_ploidy : int, optional
        Ploidy level of the genomes from the reference population. Default: 1 (phased data).
    tgt_ploidy : int, optional
        Ploidy level of the genomes from the target population. Default: 1 (phased data).
    src_ploidy : int, optional
        Ploidy level of the genomes from the source population. Default: 1 (phased data).
    out_ploidy : int, optional
        Ploidy level of the genomes from the outgroup. Default: 1 (phased data).

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        Allele frequencies for (ref, tgt, src, out) populations.
    """
    ref_freq = calc_freq(ref_gts, ref_ploidy)
    tgt_freq = calc_freq(tgt_gts, tgt_ploidy)
    src_freq = calc_freq(src_gts, src_ploidy)
    if out_gts is None:
        out_freq = np.zeros_like(ref_freq)
    else:
        out_freq = calc_freq(out_gts, out_ploidy)

    return ref_freq, tgt_freq, src_freq, out_freq


def calc_pattern_sum(
    ref_freq: np.ndarray,
    tgt_freq: np.ndarray,
    src_freq: np.ndarray,
    out_freq: np.ndarray,
    pattern: str,
) -> float:
    """
    Applies an ABBA-like pattern and returns the sum over loci of the transformed frequency products.

    Parameters
    ----------
    ref_freq:
        Allele frequencies for the reference population (no introgression) across loci.
    tgt_freq:
        Allele frequencies for the target population (receive introgression) across loci.
    src_freq:
        Allele frequencies for the source population (provide introgression) across loci.
    out_freq:
        Allele frequencies for the outgroup across loci.
    pattern : str
        A 4-character pattern string (e.g., 'abba'), where:
        - 'a': use 1 - freq
        - 'b': use freq

    Returns
    -------
    float
        Sum over loci of the product defined by the pattern.

    Raises
    ------
    ValueError
        - If the pattern string is not exactly four characters long.
        - If the pattern contains characters other than 'a' or 'b'.
    """
    if len(pattern) != 4:
        raise ValueError("Pattern must be a four-character string.")

    freqs = [ref_freq, tgt_freq, src_freq, out_freq]
    product = np.ones_like(ref_freq)

    for i, c in enumerate(pattern.lower()):
        if c == "a":
            product *= 1 - freqs[i]
        elif c == "b":
            product *= freqs[i]
        else:
            raise ValueError(
                f"Invalid character '{c}' in pattern. Only 'a' and 'b' allowed."
            )

    return float(np.sum(product))
