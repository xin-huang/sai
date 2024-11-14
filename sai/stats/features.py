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


def calc_freq(gts, ploidy=1):
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
    if ploidy == 1:
        return np.mean(gts, axis=1)
    else:
        return np.sum(gts, axis=1) / (gts.shape[1] * ploidy)


def calc_u(
    ref_gts: np.ndarray,
    tgt_gts: np.ndarray,
    src_gts_list: list[np.ndarray],
    w: float,
    x: float,
    y_list: list[float],
    ploidy: int = 1,
) -> int:
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

    Returns
    -------
    int
        The count of loci that meet all specified frequency conditions.
    """
    # Validate input parameters
    if not (0 <= w <= 1 and 0 <= x <= 1):
        raise ValueError("Parameters w and x must be within the range [0, 1].")
    if not all(0 <= y <= 1 for y in y_list):
        raise ValueError("All values in y_list must be within the range [0, 1].")
    if len(src_gts_list) != len(y_list):
        raise ValueError("The length of src_gts_list and y_list must match.")

    # Calculate allele frequencies for each group
    ref_freq = calc_freq(ref_gts, ploidy)
    tgt_freq = calc_freq(tgt_gts, ploidy)
    src_freq_list = [calc_freq(src_gts, ploidy) for src_gts in src_gts_list]

    # Check if all src_freq values match y or 1 - y for each locus
    all_match_y = np.all(
        [src_freq == y for src_freq, y in zip(src_freq_list, y_list)], axis=0
    )
    all_match_1_minus_y = np.all(
        [src_freq == 1 - y for src_freq, y in zip(src_freq_list, y_list)], axis=0
    )

    # Directly modify ref_freq and tgt_freq where src_freq matches 1 - y
    ref_freq[all_match_1_minus_y] = 1 - ref_freq[all_match_1_minus_y]
    tgt_freq[all_match_1_minus_y] = 1 - tgt_freq[all_match_1_minus_y]

    # Apply final condition: loci that match the adjusted conditions in ref and tgt
    condition = (all_match_y | all_match_1_minus_y) & (ref_freq < w) & (tgt_freq > x)

    # Count loci meeting all specified conditions
    count = np.sum(condition)
    return count


def calc_q(
    ref_gts: np.ndarray,
    tgt_gts: np.ndarray,
    src_gts_list: list[np.ndarray],
    w: float,
    y_list: list[float],
    quantile: float = 0.95,
    ploidy: int = 1,
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

    Returns
    -------
    float
        The specified quantile of the derived allele frequencies in `tgt_gts` for loci meeting the specified conditions,
        or NaN if no loci meet the criteria.
    """
    # Validate input parameters
    if not (0 <= w <= 1 and 0 <= quantile <= 1):
        raise ValueError("Parameters w and quantile must be within the range [0, 1].")
    if not all(0 <= y <= 1 for y in y_list):
        raise ValueError("All values in y_list must be within the range [0, 1].")
    if len(src_gts_list) != len(y_list):
        raise ValueError("The length of src_gts_list and y_list must match.")

    # Calculate allele frequencies for each group
    ref_freq = calc_freq(ref_gts, ploidy)
    tgt_freq = calc_freq(tgt_gts, ploidy)
    src_freq_list = [calc_freq(src_gts, ploidy) for src_gts in src_gts_list]

    # Check if all src_freq values match y or 1 - y for each locus
    all_match_y = np.all(
        [src_freq == y for src_freq, y in zip(src_freq_list, y_list)], axis=0
    )
    all_match_1_minus_y = np.all(
        [src_freq == 1 - y for src_freq, y in zip(src_freq_list, y_list)], axis=0
    )

    # Directly modify ref_freq and tgt_freq where src_freq matches 1 - y
    ref_freq[all_match_1_minus_y] = 1 - ref_freq[all_match_1_minus_y]
    tgt_freq[all_match_1_minus_y] = 1 - tgt_freq[all_match_1_minus_y]

    # Apply the final condition based on adjusted ref_freq and tgt_freq
    condition = (all_match_y | all_match_1_minus_y) & (ref_freq < w)

    # Filter `tgt_gts` frequencies based on the combined condition
    filtered_tgt_freq = tgt_freq[condition]

    # Return NaN if no loci meet the criteria
    if filtered_tgt_freq.size == 0:
        return np.nan

    # Calculate and return the specified quantile of the filtered `tgt_gts` frequencies
    return np.quantile(filtered_tgt_freq, quantile)
