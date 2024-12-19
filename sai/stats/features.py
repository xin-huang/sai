# Copyright 2024 Xin Huang and Florian R. Schmidt
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
    across reference, target, and multiple source genotypes.

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

    Raises
    ------
    ValueError
        If `w`, `x`, or any value in `y_list` is outside the range [0, 1], or if `y_list` length does not match `src_gts_list`.
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

    # Set initial condition for reference and target populations
    condition = (ref_freq < w) & (tgt_freq > x)

    # Add conditions for each source population
    for src_freq, y in zip(src_freq_list, y_list):
        condition &= src_freq == y

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
    across reference and multiple source genotypes.

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

    Raises
    ------
    ValueError
        If `w`, `quantile`, or any value in `y_list` is outside the range [0, 1], or if `y_list` length does not match `src_gts_list`.
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

    # Set initial condition for reference population
    condition = ref_freq < w

    # Add conditions for each source population
    for src_freq, y in zip(src_freq_list, y_list):
        condition &= src_freq == y

    # Filter `tgt_gts` frequencies based on the combined condition
    filtered_tgt_freq = tgt_freq[condition]

    # Return NaN if no loci meet the criteria
    if filtered_tgt_freq.size == 0:
        return np.nan

    # Calculate and return the specified quantile of the filtered `tgt_gts` frequencies
    return np.quantile(filtered_tgt_freq, quantile)


def calc_seq_div(gts1, gts2):
    """
    Calculates pairwise sequence divergence between two populations using the Hamming distance
    that supports multiallelic data (e.g., values other than 0 and 1).

    Parameters
    ----------
    gts1 : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual in the first population.
    gts2 : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual in the second population.

    Returns
    -------
    float
        The average sequence divergence between the two populations.
    """
    # Expand dimensions to broadcast `gts1` and `gts2` across each other's individuals
    expanded_gts1 = gts1[:, :, np.newaxis]  # Shape: (loci, individuals_gts1, 1)
    expanded_gts2 = gts2[:, np.newaxis, :]  # Shape: (loci, 1, individuals_gts2)

    # Calculate divergence for each pair by comparing values directly
    div_matrix = expanded_gts1 != expanded_gts2  # True where values differ

    # Average divergence across loci and individuals
    pairwise_divergence = np.sum(
        div_matrix, axis=0
    )  # Shape: (individuals_gts1, individuals_gts2)

    return pairwise_divergence


def calc_rd(ref_gts, tgt_gts, src_gts):
    """
    Calculates the average ratio of the sequence divergence between an individual from the source population
    and an individual from the admixed population, and the sequence divergence between an individual from the
    source population and an individual from the non-admixed population.

    Parameters
    ----------
    ref_gts : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual in the non-admixed population.
    tgt_gts : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual in the admixed population.
    src_gts : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual in the source population.

    Returns
    -------
    float
        The average divergence ratio.
    """
    # Step 1: Calculate sequence divergence between source and non-admixed population
    divergence_src_ref = calc_seq_div(src_gts, ref_gts)

    # Step 2: Calculate sequence divergence between source and admixed population
    divergence_src_tgt = calc_seq_div(src_gts, tgt_gts)

    if np.mean(divergence_src_ref) != 0:
        return np.mean(divergence_src_tgt) / np.mean(divergence_src_ref)
    else:
        return np.nan

    # Step 3: Replace zeros in divergence_src_ref with -1 to handle division by zero
    #divergence_src_ref_safe = np.where(
    #    divergence_src_ref == 0, np.nan, divergence_src_ref
    #)

    # Step 4: Calculate the pairwise ratios
    #divergence_ratios = (
    #    divergence_src_tgt[:, :, np.newaxis] / divergence_src_ref_safe[:, np.newaxis, :]
    #)

    # Step 5: Calculate the mean of the pairwise ratios
    #average_divergence_ratio = np.nanmean(divergence_ratios)

    #return average_divergence_ratio


def calc_abba_baba(
    ref_gts: np.ndarray,
    tgt_gts: np.ndarray,
    src_gts: np.ndarray,
    ploidy: int = 1,
) -> tuple[float, float]:
    """
    Calculates a D-statistic (ABBA-BABA) of derived allele frequencies in `tgt_gts` for loci that meet specific conditions
    across reference and multiple source genotypes.

    Parameters
    ----------
    ref_gts : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual in the reference group.
    tgt_gts : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual in the target group.
    src_gts_list : list of np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual in the target group.
    ploidy : int, optional
        The ploidy level of the organism. Default is 1, which assumes phased data.

    Returns
    -------
    tuple[float, float]
        ABBA, BABA sums of the derived allele frequencies

    """
    # Calculate allele frequence for each group
    ref_freq = calc_freq(ref_gts, ploidy)
    tgt_freq = calc_freq(tgt_gts, ploidy)
    src_freq = calc_freq(src_gts, ploidy)

    # Calculate abba/baba vectors
    abba_vec = (1.0 - tgt_freq) * src_freq * ref_freq
    baba_vec = tgt_freq * (1 - src_freq) * ref_freq
    
    # Calculate abba/baba sums
    abba = np.sum(abba_vec)
    baba = np.sum(baba_vec)

    # return ABBA-BABA sums
    return abba, baba
    

def calc_d(
        ref_gts: np.ndarray,
        tgt_gts: np.ndarray,
        src_gts: np.ndarray,
        ploidy: int = 1
) -> float:
    """
    Calculates a D-statistic (ABBA-BABA) of derived allele frequencies in `tgt_gts` for loci that meet specific conditions
    across reference and multiple source genotypes.

    Parameters
    ----------
    ref_gts : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual in the reference group.
    tgt_gts : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual in the target group.
    src_gts_list : list of np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual in the target group.
    ploidy : int, optional
        The ploidy level of the organism. Default is 1, which assumes phased data.

    Returns
    -------
    float
        D-statistic of the derived allele frequencies, or NaN if no loci meet the criteria.

    """
    # Calculate abba-baba sums with calc_abba_baba
    abba, baba = calc_abba_baba(ref_gts, tgt_gts, src_gts, ploidy)    
    
    # Return D-statistic if abba/baba greater 0
    if (abba + baba) > 0:
        return (abba - baba) / (abba + baba)
    else:
        return np.nan


def calc_fd(
        ref_gts: np.ndarray,
        tgt_gts: np.ndarray,
        src_gts: np.ndarray,
        ploidy: int = 1
):
    """
    Calculates a fD-statistic of derived allele frequencies in `tgt_gts` for loci that meet specific conditions
    across reference and multiple source genotypes.

    Parameters
    ----------
    ref_gts : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual in the reference group.
    tgt_gts : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual in the target group.
    src_gts_list : list of np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual in the target group.
    ploidy : int, optional
        The ploidy level of the organism. Default is 1, which assumes phased data.

    Returns
    -------
    float
        D-statistic of the derived allele frequencies, or NaN if no loci meet the criteria.

    """
    # Calculate allele frequence for each group
    ref_freq = calc_freq(ref_gts, ploidy)
    tgt_freq = calc_freq(tgt_gts, ploidy)
    src_freq = calc_freq(src_gts, ploidy)

    # Calculate abba/baba fd vectors
    abba_fd1_vec = (1.0 - tgt_freq) * src_freq * src_freq
    baba_fd1_vec = tgt_freq * (1 - src_freq) * src_freq
    abba_fd2_vec = (1.0 - tgt_freq) * ref_freq * ref_freq
    baba_fd2_vec = tgt_freq * (1 - ref_freq) * ref_freq
    
    # ref and src conditions
    fd1_condition = (src_freq > ref_freq)
    fd2_condition = (src_freq < ref_freq)
    
    # Calculate fd abba/baba sums
    abba_fd = fd1_condition * abba_fd1_vec + fd2_condition * abba_fd2_vec
    abba_fd = np.sum(abba_fd)
    
    baba_fd = fd1_condition * baba_fd1_vec + fd2_condition * baba_fd2_vec
    baba_fd = np.sum(baba_fd)
    
    # Return fD-statistic if abba/baba greater 0
    if (abba_fd + baba_fd) > 0:
        # Calculate normal abba, baba sum values from calc_abba_baba
        abba, baba = calc_abba_baba(ref_gts, tgt_gts, src_gts, ploidy)
        
        return (abba - baba) / (abba_fd + baba_fd)
    else:
        return np.nan
    