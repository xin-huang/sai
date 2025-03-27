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
    """
    return np.sum(gts, axis=1) / (gts.shape[1] * ploidy)


def compute_matching_loci(
    ref_gts: np.ndarray,
    tgt_gts: np.ndarray,
    src_gts_list: list[np.ndarray],
    w: float,
    y_list: list[tuple[str, float]],
    ploidy: int,
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
    ploidy : int
        The ploidy level of the organism.
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
    ref_freq = calc_freq(ref_gts, ploidy)
    tgt_freq = calc_freq(tgt_gts, ploidy)
    src_freq_list = [calc_freq(src_gts, ploidy) for src_gts in src_gts_list]

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


def calc_u(
    ref_gts: np.ndarray,
    tgt_gts: np.ndarray,
    src_gts_list: list[np.ndarray],
    pos: np.ndarray,
    w: float,
    x: float,
    y_list: list[float],
    ploidy: int = 1,
    anc_allele_available: bool = False,
) -> tuple[int, np.ndarray]:
    """
    Calculates the count of genetic loci that meet specified allele frequency conditions
    across reference, target, and multiple source genotypes, with adjustments based on src_freq consistency.

    Parameters
    ----------
    ref_gts : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual in the reference group.
    tgt_gts : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual in the target group.
    src_gts_list : list of np.ndarray
        A list of 2D numpy arrays for each source population, where each row represents a locus and each column
        represents an individual in that source population.
    pos : np.ndarray
        A 1D numpy array where each element represents the genomic position.
    w : float
        Threshold for the allele frequency in `ref_gts`. Only loci with frequencies less than `w` are counted.
        Must be within the range [0, 1].
    x : float
        Threshold for the allele frequency in `tgt_gts`. Only loci with frequencies greater than `x` are counted.
        Must be within the range [0, 1].
    y_list : list of float
        List of exact allele frequency thresholds for each source population in `src_gts_list`.
        Must be within the range [0, 1] and have the same length as `src_gts_list`.
    ploidy : int, optional
        The ploidy level of the organism. Default is 1, which assumes phased data.
    anc_allele_available : bool
        If True, checks only for matches with `y` (assuming `1` represents the derived allele).
        If False, checks both matches with `y` and `1 - y`, taking the major allele in the source as the reference.

    Returns
    -------
    tuple[int, np.ndarray]
        - The count of loci that meet all specified frequency conditions.
        - A 1D numpy array containing the genomic positions of the loci that meet the conditions.

    Raises
    ------
    ValueError
        If `x` is outside the range [0, 1].
    """
    # Validate input parameters
    if not (0 <= x <= 1):
        raise ValueError("Parameter x must be within the range [0, 1].")

    ref_freq, tgt_freq, condition = compute_matching_loci(
        ref_gts,
        tgt_gts,
        src_gts_list,
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
    return count, loci_positions


def calc_q(
    ref_gts: np.ndarray,
    tgt_gts: np.ndarray,
    src_gts_list: list[np.ndarray],
    pos: np.ndarray,
    w: float,
    y_list: list[float],
    quantile: float = 0.95,
    ploidy: int = 1,
    anc_allele_available: bool = False,
) -> float:
    """
    Calculates a specified quantile of derived allele frequencies in `tgt_gts` for loci that meet specific conditions
    across reference and multiple source genotypes, with adjustments based on src_freq consistency.

    Parameters
    ----------
    ref_gts : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual in the reference group.
    tgt_gts : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual in the target group.
    src_gts_list : list of np.ndarray
        A list of 2D numpy arrays for each source population, where each row represents a locus and each column
        represents an individual in that source population.
    pos: np.ndarray
        A 1D numpy array where each element represents the genomic position.
    w : float
        Frequency threshold for the derived allele in `ref_gts`. Only loci with frequencies lower than `w` are included.
        Must be within the range [0, 1].
    y_list : list of float
        List of exact frequency thresholds for each source population in `src_gts_list`.
        Must be within the range [0, 1] and have the same length as `src_gts_list`.
    quantile : float, optional
        The quantile to compute for the filtered `tgt_gts` frequencies. Must be within the range [0, 1].
        Default is 0.95 (95% quantile).
    ploidy : int, optional
        The ploidy level of the organism. Default is 1, which assumes phased data.
    anc_allele_available : bool
        If True, checks only for matches with `y` (assuming `1` represents the derived allele).
        If False, checks both matches with `y` and `1 - y`, taking the major allele in the source as the reference.

    Returns
    -------
    tuple[float, np.ndarray]
        - The specified quantile of the derived allele frequencies in `tgt_gts` for loci meeting the specified conditions,
          or NaN if no loci meet the criteria.
        - A 1D numpy array containing the genomic positions of the loci that meet the conditions.

    Raises
    ------
    ValueError
        If `quantile` is outside the range [0, 1].
    """
    # Validate input parameters
    if not (0 <= quantile <= 1):
        raise ValueError("Parameter quantile must be within the range [0, 1].")

    ref_freq, tgt_freq, condition = compute_matching_loci(
        ref_gts,
        tgt_gts,
        src_gts_list,
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
        return np.nan, np.array([])

    threshold = np.nanquantile(filtered_tgt_freq, quantile)
    loci_positions = filtered_positions[filtered_tgt_freq >= threshold]

    # Calculate and return the specified quantile of the filtered `tgt_gts` frequencies
    return threshold, loci_positions
