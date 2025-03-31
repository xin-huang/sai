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
from typing import Tuple, Union
from scipy.spatial.distance import cdist
from itertools import product


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
    if ploidy == 1:
        return np.mean(gts, axis=1)
    else:
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


def calc_rd(
    src_gts: np.ndarray,
    tgt_gts: np.ndarray,
    ref_gts: np.ndarray = None,
    metric: str = "hamming",
) -> float:
    """
    Compute the average ratio of sequence divergence between an individual from the
    source population and an individual from the target (admixed) population, compared
    to the sequence divergence between an individual from the source population and
    an individual from the reference (non-admixed) population.

    This is done by comparing the pairwise distances between individuals in the source
    group (`src_gts`) with those in the target group (`tgt_gts`) and the reference group
    (`ref_gts`) using a chosen distance metric. The mean of the pairwise distances for
    each row in both matrices is computed, and the ratio of these means is averaged over
    all combinations of source-target and source-reference rows.

    Parameters
    ----------
    src_gts : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents
        an individual in the source population.

    tgt_gts : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents
        an individual in the target (admixed) population.

    ref_gts : np.ndarray, optional
        A 2D numpy array where each row represents a locus and each column represents
        an individual in the reference (non-admixed) population. If `None`, a zero array
        of the same shape as `tgt_gts` is used. Default is `None`.

    metric : str, optional
        The distance metric to use for the pairwise distance calculation. Default is "hamming".

    Returns
    -------
    float
        The computed average ratio of sequence divergence between the source-target and
        source-reference populations.

    """

    # If ref_gts is None, set it to a zero matrix of the same shape as tgt_gts
    if ref_gts is None:
        ref_gts = np.zeros(tgt_gts.shape)

    # pairwise distances
    seq_divs_src_tgt = cdist(src_gts.T, tgt_gts.T, metric=metric)
    seq_divs_src_ref = cdist(src_gts.T, ref_gts.T, metric=metric)

    # mean of each row
    mean_src_tgt = np.mean(seq_divs_src_tgt, axis=1)
    mean_src_ref = np.mean(seq_divs_src_ref, axis=1)

    all_pairs = list(
        product(range(seq_divs_src_tgt.shape[0]), range(seq_divs_src_ref.shape[0]))
    )
    count = len(all_pairs)

    average_r = 0

    # Loop over all combinations
    for i, j in all_pairs:
        # Use precomputed row means
        row_tgt_mean = mean_src_tgt[i]
        row_ref_mean = mean_src_ref[j]

        # Accumulate the ratio of means
        if row_ref_mean != 0:
            average_r += row_tgt_mean / row_ref_mean
        else:
            print("Warning! An average ratio is zero!")

    average_r = average_r / count
    return average_r


# functions taken from LogisticRegression/gaia


def cal_n_ton(tgt_gt, is_phased, ploidy):
    """
    Description:
        Calculates individual frequency spectra for samples.

    Arguments:
        tgt_gt numpy.ndarray: Genotype matrix from the target population.
        ploidy int: Ploidy of the genomes.

    Returns:
        spectra numpy.ndarray: Individual frequency spectra for haplotypes.
    """
    if is_phased:
        ploidy = 1
    mut_num, sample_num = tgt_gt.shape
    iv = np.ones((sample_num, 1))
    counts = (tgt_gt > 0) * np.matmul(tgt_gt, iv)
    spectra = np.array(
        [
            np.bincount(
                counts[:, idx].astype("int64"), minlength=sample_num * ploidy + 1
            )
            for idx in range(sample_num)
        ]
    )
    # ArchIE does not count non-segragating sites
    spectra[:, 0] = 0

    return spectra


def cal_dist(gt1, gt2):
    """
    Description:
        Calculates pairwise Euclidean distances between two genotype matrixes.

    Arguments:
        gt1 numpy.ndarray: Genotype matrix 1.
        gt2 numpy.ndarray: Genotype matrix 2.

    Returns:
        dists numpy.ndarray: Distances estimated.
    """

    from scipy.spatial import distance_matrix

    dists = distance_matrix(np.transpose(gt2), np.transpose(gt1))
    dists.sort()

    return dists


# helper features


def return_segsite_occurrence_vector(
    tgt_gt, ploidy: int = 1, remove_non_segregating: bool = False
):
    """
    Compute the number of segregating sites where the mutant allele occurs exactly i times in the sample.

    Arguments:
        tgt_gt (numpy.ndarray): Genotype matrix (mutations x samples).
        remove_non_segregating (bool): if True, return a vector which only proviedes the counts for at least one segregating site, i.e. the first entry is removed

    Returns:
        seg_sites (numpy.ndarray): A vector where index i represents the count of sites where
                                   the mutant allele appears exactly i times.

    """
    # Compute the sum of mutant alleles per mutation (sum across rows)
    mutant_counts = np.sum(tgt_gt, axis=1)

    seg_sites = np.bincount(mutant_counts, minlength=tgt_gt.shape[1] * ploidy)
    if remove_non_segregating:
        seg_sites = seg_sites[1:]

    return seg_sites


# MaLAdapt features


def heterozygosity(tgt_gts: np.ndarray, ploidy: int = 1) -> float:
    """
    Computes the expected heterozygosity (H) for a given genotype matrix.

    Expected heterozygosity is calculated as:
        H = mean(2 * p * (1 - p)),
    where p is the allele frequency at each locus.

    Parameters
    ----------
    tgt_gts : np.ndarray
        A 2D numpy array where each row represents a locus and each column represents an individual.
    ploidy : int, optional
        The ploidy level of the organism (default is 1).

    Returns
    -------
    float
        The mean expected heterozygosity across all loci.
    """
    gts_freq = calc_freq(tgt_gts, ploidy=ploidy)
    hetvec = 2 * gts_freq * (1.0 - gts_freq)
    Het = np.mean(hetvec)
    return Het


def num_segregating_sites(
    gts: np.ndarray, return_frequencies: bool = False
) -> Union[int, tuple[int, np.ndarray]]:
    """
    Computes the number of segregating sites in a genotype matrix.

    A site (locus) is considered segregating if it contains at least two different alleles (i.e.,
    the allele frequency is between 0 and 1).

    Parameters
    ----------
    gts : np.ndarray
        A 2D NumPy array where each row represents a locus and each column represents an individual.
    return_frequencies : bool, optional
        If True, also returns the allele frequencies at segregating sites. Default is False.

    Returns
    -------
    int
        The number of segregating sites.
    tuple[int, np.ndarray], optional
        If `return_frequencies` is True, returns a tuple containing the number of segregating sites
        and an array of their frequencies.
    """
    # Compute allele frequencies for each locus
    gts_freq = calc_freq(gts)

    # Identify segregating sites (where 0 < p < 1)
    segregating_mask = (gts_freq > 0) & (gts_freq < 1)

    num_S = np.sum(segregating_mask)

    if return_frequencies:
        return num_S, gts_freq[segregating_mask]
    return num_S


def S_ABBA_BABA_calc(
    pop1_freq: np.ndarray,
    pop2_freq: np.ndarray,
    pop3_freq: np.ndarray,
    outgroup_freq: float = 0,
    return_numerator_for_S: bool = False,
) -> float:
    """
    Calculate the ABBA-BABA statistic, which is used in population genetics to test
    for asymmetry in allele sharing between populations and an outgroup. The function
    computes either the D-statistic (ABBA - BABA) / (ABBA + BABA) or the numerator
    (ABBA - BABA) based on the specified parameters.

    Parameters
    ----------
    pop1_freq : np.ndarray
        An array of allele frequencies for population 1 at each locus. (Usually the reference population)

    pop2_freq : np.ndarray
        An array of allele frequencies for population 2 at each locus. (Usually the target population)

    pop3_freq : np.ndarray
        An array of allele frequencies for population 3 at each locus. (Usually the source population)

    outgroup_freq : float, optional, default=0
        The allele frequency in the outgroup population. The default is 0.

    return_numerator_for_S : bool, optional, default=False
        If True, the function returns the numerator (ABBA - BABA) instead of
        the D-statistic. The default is False.

    Returns
    -------
    float
        The D-statistic (if `return_numerator_for_S` is False) or the numerator
        (ABBA - BABA) if `return_numerator_for_S` is True. Returns `nan` if the
        denominator is zero.

    Notes
    -----
    The ABBA-BABA statistic is commonly used in studies of introgression or
    gene flow between populations, where a high value of the D-statistic suggests
    gene flow from population 1 to 2, and a low value suggests gene flow in the
    opposite direction.
    """

    abbavec = (1.0 - pop1_freq) * pop2_freq * pop3_freq * (1 - outgroup_freq)
    babavec = pop1_freq * (1.0 - pop2_freq) * pop3_freq * (1 - outgroup_freq)

    # Summing up across loci
    abba = np.sum(abbavec)
    baba = np.sum(babavec)

    if not return_numerator_for_S:
        # Compute D-statistic
        if (abba + baba) > 0:
            return (abba - baba) / (abba + baba)
        else:
            return float("nan")
    else:
        return abba - baba


def compute_ABBA_BABA_D(
    src_gts: np.ndarray,
    ref_gts: np.ndarray,
    tgt_gts: np.ndarray,
    out_gts: np.ndarray = None,
    ploidy: int = 1,
) -> float:
    """
    Computes Patterson's D-statistic (ABBA-BABA statistic) for detecting admixture between populations.

    Parameters
    ----------
    src_gts : np.ndarray
        A 2D numpy array representing haplotypes of population 1 (source).
    ref_gts : np.ndarray
        A 2D numpy array representing haplotypes of population 2 (reference / sister group).
    tgt_gts : np.ndarray
        A 2D numpy array representing haplotypes of population 3 (target).
    out_gts : np.ndarray
        A 2D numpy array representing haplotypes of population 4 (outgroup).
        If not provided, it is assumed that the ancestral allel is always present in the outgroup, and thus the frequency of the derived allel (p4_freq) is 0.

    Returns
    -------
    float
        The D-statistic value, indicating the degree of allele sharing bias.
    """

    # Compute allele frequencies using the provided calc_freq function
    src_freq = calc_freq(src_gts, ploidy=ploidy)
    ref_freq = calc_freq(ref_gts, ploidy=ploidy)
    tgt_freq = calc_freq(tgt_gts, ploidy=ploidy)
    if out_gts is not None:
        out_freq = calc_freq(out_gts, ploidy=ploidy)
    else:
        out_freq = 0

    D_stat = S_ABBA_BABA_calc(
        ref_freq,
        tgt_freq,
        src_freq,
        outgroup_freq=out_freq,
        return_numerator_for_S=False,
    )

    return D_stat


def compute_fd(
    src_gts: np.ndarray,
    ref_gts: np.ndarray,
    tgt_gts: np.ndarray,
    out_gts: np.ndarray = None,
    ploidy: int = 1,
    use_hom: bool = False,
) -> float:
    """
    Computes fD-statistic for detecting admixture between populations (Martin 2015).

    Parameters
    ----------

    src_gts : np.ndarray
        A 2D numpy array representing haplotypes of population 3 (source).
    ref_gts : np.ndarray
        A 2D numpy array representing haplotypes of population 1 (reference / sister group).
    tgt_gts : np.ndarray
        A 2D numpy array representing haplotypes of population 2 (target).
    out_gts : np.ndarray
        A 2D numpy array representing haplotypes of population 4 (outgroup).
        If not provided, it is assumed that the ancestral allel is always present in the outgroup, and thus the frequency of the derived allel (p4_freq) is 0.
    use_hom: boolean
        If true, compute fhom instead of fd (Martin 2015, p.254)

    Returns
    -------
    float
        The fd-statistic value.
    """

    # Compute allele frequencies using the provided calc_freq function
    src_freq = calc_freq(src_gts)  # 3
    ref_freq = calc_freq(ref_gts)  # 1
    tgt_freq = calc_freq(tgt_gts)  # 2
    if out_gts is not None:
        out_freq = calc_freq(out_gts, ploidy=ploidy)
    else:
        out_freq = 0

    fd_nominator = S_ABBA_BABA_calc(
        ref_freq,
        tgt_freq,
        src_freq,
        outgroup_freq=out_freq,
        return_numerator_for_S=True,
    )

    if not use_hom:
        donor_pop_freq = np.maximum(src_freq, tgt_freq)
        fd_denumerator = S_ABBA_BABA_calc(
            ref_freq,
            donor_pop_freq,
            donor_pop_freq,
            outgroup_freq=out_freq,
            return_numerator_for_S=True,
        )
    else:
        fd_denumerator = S_ABBA_BABA_calc(
            ref_freq,
            tgt_freq,
            tgt_freq,
            outgroup_freq=out_freq,
            return_numerator_for_S=True,
        )

    fd = fd_nominator / fd_denumerator
    return fd


def theta_W(gts: np.ndarray, ploidy: int = 1) -> float:
    """
    Calculates Watterson's Theta from genotype data.

    Parameters
    ----------
    gts : np.ndarray
        A 2D numpy array where rows represent genetic sites and columns represent individuals.
        Each element in the array represents the allele at a given site for an individual, typically
        0 (reference allele) or 1 (alternate allele).

    Returns
    -------
    float
        theta_W
    """
    S = num_segregating_sites(gts)
    individuals = gts.shape[1] * ploidy
    harmonic_number = sum(1 / i for i in range(1, individuals + 1))
    theta_W = S / harmonic_number
    return theta_W


def theta_pi_maladapt(gts: np.ndarray, ploidy: int = 1) -> float:
    """
    Calculates theta_pi from genotype data as in MaLAdapt.

    Parameters
    ----------
    gts : np.ndarray
        A 2D numpy array where rows represent genetic sites and columns represent individuals.
    ploidy : int, optional
        The number of chromosome copies per individual (default is 1 for haploid).

    Returns
    -------
    float
        theta_pi, estimator  of nucleotide diversity.
    """
    # Compute segregating sites (S) and allele frequencies (S_freq)
    S, S_freq = num_segregating_sites(gts, return_frequencies=True)

    # Total number of alleles
    individuals = gts.shape[1] * ploidy

    # Compute π using allele frequencies
    pi = sum(2 * S_freq * (1.0 - S_freq))

    # Adjust for sample size
    theta_pi_value = pi * individuals / (individuals - 1)

    return theta_pi_value


def theta_pi(gts: np.ndarray, ploidy: int = 1) -> float:
    """
    Compute theta_pi from a genotype matrix,
    following Zeng et al. 2006: Statistical Tests for Detecting Positive Selection by Utilizing High-Frequency Variants

    The function calculates θπ (nucleotide diversity) using a combination of the site frequency
    spectrum (SFS) and pairwise comparisons between individuals. It normalizes the result by the
    total number of possible pairwise comparisons, given the ploidy of the individuals.
    Should in principle also work for unphased data.

    Parameters
    ----------
    gts : np.ndarray
        A 2D NumPy array of shape (mutations x samples) where each entry represents the number of
        mutant alleles (0, 1, or 2, ...) at each site for each individual.

    ploidy : int, optional, default=1
        The ploidy of the population. For example, ploidy = 2 for diploid organisms and ploidy = 1 for haploid organisms.

    Returns
    -------
    float
        The nucleotide diversity (θπ) for the given genotype matrix, considering the number of pairwise comparisons.

    """
    from math import comb

    individuals = gts.shape[1] * ploidy
    xi = return_segsite_occurrence_vector(
        gts, remove_non_segregating=False, ploidy=ploidy
    )

    combination_prefactor = 1 / (comb(individuals, 2))

    pi_sum = 0
    # the first term (all ancestral) is automatically 0, the last term (all derived) we remove from the calculation
    for ie, entry in enumerate(xi[:-1]):
        pi_sum = pi_sum + ie * (individuals - ie) * entry

    theta_pi = pi_sum * combination_prefactor

    return theta_pi


def compute_Fay_Wu_theta_h(gts: np.ndarray, ploidy: int = 1) -> float:
    """
    Computes Fay and Wu's theta_H based on genotype matrix (gts).

    This function calculates theta_H, a measure of nucleotide diversity that
    accounts for the number of segregating sites in a given genotype matrix.
    Theta_H is commonly used in population genetics to estimate genetic diversity.

    Parameters
    ----------
    gts : np.ndarray
        A 2D numpy array where rows represent loci and columns represent individuals.
        The matrix contains genotype information, where alleles are coded numerically.

    ploidy : int, optional, default=1
        The ploidy level of the individuals. For diploid organisms, ploidy would be 2.

    Returns
    -------
    float
        The calculated Fay and Wu's theta_H value.

    Notes
    -----
    Theta_H is calculated using the formula:
    theta_H = (2 * sum(x_vals^2 * xi[:-1])) / (n * (n - 1))
    where:
        - xi is the vector of derived allele counts,
        - x_vals is the range of indices for segregating sites,
        - n is the total haploid sample size.

    """

    individuals = gts.shape[1] * ploidy
    xi = return_segsite_occurrence_vector(
        gts, remove_non_segregating=False, ploidy=ploidy
    )

    x_vals = np.arange(len(xi) - 1)

    num = 2 * x_vals**2 * xi[:-1]
    denom = individuals * (individuals - 1)

    theta_h = np.sum(num) / denom
    return theta_h


def H1_H12_values(
    gts: np.ndarray,
    only_derived_homozygous: bool = False,
    compute_H123: bool = False,
    ploidy: int = 1,
) -> Union[Tuple[float, float], Tuple[float, float, float]]:
    """
    Computes H1, H12, and optionally H123 statistics to measure haplotype homozygosity.
    Garud 2015

    Parameters
    ----------
    gts : np.ndarray
        Genotype data (allele frequencies or haplotype counts).
    only_derived_homozygous : bool, optional
        If True, considers only derived homozygous frequencies (default is False).
    compute_H123 : bool, optional
        If True, also computes the H123 statistic (default is False).
    ploidy : int, optional
        Number of copies of each chromosome per individual (default is 1).

    Returns
    -------
    Tuple[float, float] or Tuple[float, float, float]
        H1 and H12 values; if compute_H123 is True, returns H1, H12, and H123 values.
    """
    freqs = calc_freq(gts, ploidy=ploidy)
    if not only_derived_homozygous:
        freqs[freqs < 0.5] = 1 - freqs[freqs < 0.5]
    freqs = np.sort(freqs, axis=0)[::-1]

    H1_value = np.sum(freqs**2)

    # alternative calculation
    # H12_value_alt = H1_value + 2 * freqs[0] * freqs[1]

    H12_value = ((freqs[0] + freqs[1]) ** 2) + np.sum(freqs[2:] ** 2)

    H123_value = ((freqs[0] + freqs[1] + freqs[2]) ** 2) + np.sum(freqs[3:] ** 2)

    if not compute_H123:
        return H1_value, H12_value
    else:
        return H1_value, H12_value, H123_value


def H2_value(
    gts: np.ndarray, only_derived_homozygous: bool = False, ploidy: int = 1
) -> float:
    """
    Computes the H2 statistic, which measures haplotype homozygosity excluding the most common haplotype.
    Garud 2015

    Parameters
    ----------
    gts : np.ndarray
        Genotype data (allele frequencies or haplotype counts).
    only_derived_homozygous : bool, optional
        If True, considers only derived homozygous frequencies (default is False).
    ploidy : int, optional
        Number of copies of each chromosome per individual (default is 1).

    Returns
    -------
    float
        The H2 value.
    """
    freqs = calc_freq(gts, ploidy=ploidy)
    if not only_derived_homozygous:
        freqs[freqs < 0.5] = 1 - freqs[freqs < 0.5]
    freqs = np.sort(freqs, axis=0)[::-1]

    H2 = ((freqs[1]) ** 2) + np.sum(freqs[2:] ** 2)
    return H2


def H2_H1_ratio(
    gts: np.ndarray, only_derived_homozygous: bool = False, ploidy: int = 1
) -> float:
    """
    Computes the H2/H1 ratio, a measure of haplotype diversity relative to the most common haplotype.
    Garud 2015

    Parameters
    ----------
    gts : np.ndarray
        Genotype data (allele frequencies or haplotype counts).
    only_derived_homozygous : bool, optional
        If True, considers only derived homozygous frequencies (default is False).
    ploidy : int, optional
        Number of copies of each chromosome per individual (default is 1).

    Returns
    -------
    float
        The H2/H1 ratio.
    """
    H2 = H2_value(gts, only_derived_homozygous=only_derived_homozygous, ploidy=ploidy)

    H1, H12 = H1_H12_values(
        gts, only_derived_homozygous=only_derived_homozygous, ploidy=1
    )

    H2_H1 = H2 / H1
    return H2_H1


def compute_LD_D_estimate(gts, ploidy=1, phased=True, compute_r=True):
    """
    Computes the linkage disequilibrium coefficient D for all pairs of SNPs in phased data.
    Using the formulae from Ragsdale 2019: Unbiased Estimation of Linkage Disequilibrium from Unphased
    Data

    Parameters:
    gts (numpy.ndarray): A 2D numpy array where rows represent SNPs and columns represent haplotypes.


    Returns:
    numpy.ndarray: A symmetric matrix where entry (i, j) is the D value between SNP i and SNP j.
    """

    if phased or ploidy == 1:
        num_snps = gts.shape[0]

        # Initialize LD matrix with NaN
        ld_matrix = np.zeros((num_snps, num_snps))
        r_matrix = np.zeros((num_snps, num_snps))

        num_snps, num_individuals = gts.shape

        # Initialize LD matrix with NaN
        ld_matrix = np.full((num_snps, num_snps), np.nan)
        if compute_r:
            r_matrix = np.zeros((num_snps, num_snps))

        for i in range(num_snps):
            for j in range(i + 1, num_snps):  # Compute only upper triangle (symmetry)
                geno1 = gts[i, :]
                geno2 = gts[j, :]

                g1, g2, g3, g4 = 0, 0, 0, 0

                for ind in range(num_individuals):
                    loci1 = geno1[ind]
                    loci2 = geno2[ind]

                    if loci1 == 0 and loci2 == 0:
                        g4 += 1
                    elif loci1 == 0 and loci2 == 1:
                        g3 += 1
                    elif loci1 == 1 and loci2 == 0:
                        g2 += 1
                    elif loci1 == 1 and loci2 == 1:
                        g1 += 1

                num1 = g1 * g4
                num2 = g2 * g3

                # -1 used as in Ragsdale 2019
                denom = num_individuals * (num_individuals - 1)

                D_comp = (num1 - num2) / denom

                ld_matrix[i, j] = D_comp
                ld_matrix[j, i] = D_comp

                p = (g1 + g2) / num_individuals
                q = (g1 + g3) / num_individuals

                if compute_r:
                    r_value = D_comp / (np.sqrt(p * (1 - p) * q * (1 - q)))
                    r_matrix[i, j] = r_value
                    r_matrix[j, i] = r_value

        np.fill_diagonal(ld_matrix, np.nan)

    else:
        raise Exception("This function is only appropriate for phased/haploid data!")
    if not compute_r:
        return ld_matrix
    else:
        return ld_matrix, r_matrix


# LD specific


def compute_LD_D(
    gts: np.ndarray,
    ploidy: int = 1,
    filter_unique: bool = True,
    compute_r: bool = True,
    phased: bool = True,
    maladapt_correction: bool = False,
) -> np.ndarray:
    """
    Computes the linkage disequilibrium coefficient D for all pairs of SNPs in phased or unphased data.

    This function computes the coefficient D for pairs of SNPs either using phased or unphased data.

    Parameters
    ----------
    gts : numpy.ndarray
        A 2D numpy array where rows represent SNPs and columns represent haplotypes (genotypes).
        The array should have shape (n_snps, n_haplotypes), where `n_snps` is the number of SNPs and
        `n_haplotypes` is the number of haplotypes (usually 2 * n_individuals for diploid data).

    ploidy : int, optional, default=1
        The ploidy of the species (e.g., 1 for haploid or 2 for diploid). This is relevant for determining
        the method of computing linkage disequilibrium (LD).

    filter_unique : bool, optional, default=True
        If True, the function will filter out rows of `gts` that contain only unique values (i.e., no variation).

    compute_r : bool, optional, default=True
        If True, the function will also compute the r coefficient alongside the D coefficient. Otherwise,
        only the D coefficient is computed.

    phased : bool, optional, default=True
        If True, the data is considered to be phased. Otherwise, it is unphased.

    maladapt_correction : bool, optional, default=False
        If True, the function will calculate LD very similar to the implementation in MaLAdapt.

    Returns
    -------
    numpy.ndarray
        A symmetric matrix where entry (i, j) is the D value between SNP i and SNP j. The matrix will have
        shape (n_snps, n_snps), where `n_snps` is the number of SNPs in the `gts` array.

    Notes
    -----
    - The function dispatches to different methods based on the ploidy and phased status of the data.
    - The function assumes the input genotypic data (`gts`) is properly formatted and may raise errors if the
      data structure does not conform to expected shapes or values.
    """

    if filter_unique:
        gts = gts[np.array([len(np.unique(row)) > 1 for row in gts])]

    if (phased or ploidy == 1) and not maladapt_correction:
        return compute_LD_D_haploid(gts, ploidy=1, compute_r=compute_r, phased=True)

    elif ploidy == 2:
        # for unphased data
        return compute_ld_burrows(gts, compute_r=compute_r)
    else:
        return compute_LD_D_maladapt(
            gts, ploidy=1, compute_r=compute_r, maladapt_correction=maladapt_correction
        )


def compute_LD_D_haploid(gts, ploidy=1, compute_r=True, phased=True):
    """
    Computes the linkage disequilibrium coefficient D for all pairs of SNPs in phased/haploid data.

    Parameters:
    gts (numpy.ndarray): A 2D numpy array where rows represent SNPs and columns represent haplotypes.


    Returns:
    numpy.ndarray: A symmetric matrix where entry (i, j) is the D value between SNP i and SNP j.
    """

    if phased or ploidy == 1:
        num_snps = gts.shape[0]

        # Initialize LD matrix with NaN
        ld_matrix = np.zeros((num_snps, num_snps))
        r_matrix = np.zeros((num_snps, num_snps))

        for i in range(num_snps):
            for j in range(i + 1, num_snps):  # Compute only upper triangle (symmetry)
                hap1 = gts[i, :]
                hap2 = gts[j, :]

                # Compute allele frequencies
                pA = np.mean(hap1)  # Frequency of major allele at SNP i
                pB = np.mean(hap2)  # Frequency of major allele at SNP j

                # Compute haplotype frequency P_AB
                pAB = np.mean(hap1 * hap2)  # Joint probability

                # Compute LD coefficient D
                D = pAB - (pA * pB)

                ld_matrix[i, j] = D
                ld_matrix[j, i] = D

                if compute_r:
                    r_value = D / (np.sqrt(pA * (1 - pA) * pB * (1 - pB)))
                    r_matrix[i, j] = r_value
                    r_matrix[j, i] = r_value

        # Optionally set diagonal to 0
        np.fill_diagonal(ld_matrix, np.nan)

    else:
        raise Exception("This function works only with phased/haploid data correctly!")

    if not compute_r:
        return ld_matrix
    else:
        return ld_matrix, r_matrix


def compute_ld_burrows(
    gts: np.ndarray,
    compute_r: bool = True,
    ploidy: int = 1,
) -> Union[Tuple[np.ndarray, np.ndarray], np.ndarray]:
    """
    Computes the linkage disequilibrium coefficient D for all pairs of SNPs in phased data.
    Using the formulae for D estimation from Ragsdale 2019: Unbiased Estimation of Linkage Disequilibrium from Unphased
    Data (Burrows “composite covariance measure of LD)
    Should work ONLY for diploid unphased!!!

    Parameters
    ----------
    gts : numpy.ndarray
        A 2D numpy array where each row represents a SNP, and each column represents a genotype for an individual.
        The shape of the array is (n_snps, n_individuals), where `n_snps` is the number of SNPs and
        `n_individuals` is the number of individuals or haplotypes.

    ploidy : int, optional, default=1
    The ploidy level, indicating the number of chromosome sets.

    compute_r : bool, optional, default=True
        If True, the function will compute the r correlation coefficient in addition to the D coefficient.
        If False, only the D coefficient will be computed.

    Returns
    -------
    tuple of numpy.ndarray or numpy.ndarray
        - If `compute_r` is True, a tuple of two numpy arrays is returned:
        - `ld_matrix`: A symmetric matrix where entry (i, j) is the D value between SNP i and SNP j.
        - `r_matrix`: A symmetric matrix where entry (i, j) is the r value between SNP i and SNP j, or NaN if the r value
            cannot be computed.
        - If `compute_r` is False, only the `ld_matrix` is returned, which contains the D values.

    """

    num_snps, num_individuals = gts.shape

    # Initialize LD matrix with NaN
    ld_matrix = np.full((num_snps, num_snps), np.nan)
    if compute_r:
        r_matrix = np.zeros((num_snps, num_snps))

    for i in range(num_snps):
        for j in range(i + 1, num_snps):  # Compute only upper triangle (symmetry)
            geno1 = gts[i, :]
            geno2 = gts[j, :]

            g1, g2, g3, g4, g5, g6, g7, g8, g9 = 0, 0, 0, 0, 0, 0, 0, 0, 0

            for ind in range(num_individuals):
                loci1 = geno1[ind]
                loci2 = geno2[ind]

                if loci1 == 0 and loci2 == 0:
                    g9 += 1
                elif loci1 == 0 and loci2 < ploidy:
                    g8 += 1
                elif loci1 == 0 and loci2 == ploidy:
                    g7 += 1

                elif loci1 < ploidy and loci2 == 0:
                    g6 += 1
                elif loci1 < ploidy and loci2 < ploidy:
                    g5 += 1
                elif loci1 < ploidy and loci2 == ploidy:
                    g4 += 1

                elif loci1 == ploidy and loci2 == 0:
                    g3 += 1
                elif loci1 == ploidy and loci2 < ploidy:
                    g2 += 1
                elif loci1 == ploidy and loci2 == ploidy:
                    g1 += 1

            g1 = g1 / num_individuals
            g2 = g2 / num_individuals
            g3 = g3 / num_individuals
            g4 = g4 / num_individuals
            g5 = g5 / num_individuals
            g6 = g6 / num_individuals
            g7 = g7 / num_individuals
            g8 = g8 / num_individuals
            g9 = g9 / num_individuals

            # alternative computation
            """       
            xAB = g1 + (g2/2) + (g4/2) + (g5/4)
            xab = g9 + (g8/2) + (g6/2) + (g5/4)
            xAb = g3 + (g2/2) + (g5/4) + (g6/2)
            xaB = g7 + (g4/2) + (g5/4) + (g8/2)
            """

            p = (g1 + g2 + g3) + ((1 / 2) * (g4 + g5 + g6))
            q = (g1 + g4 + g7) + ((1 / 2) * (g2 + g5 + g8))

            D_comp = (2 * g1 + g2 + g4 + (1 / 2) * g5) - 2 * p * q
            # alternative computation
            # D_comp = 2 * (xAB*xab-xAb*xaB)

            ld_matrix[i, j] = D_comp
            ld_matrix[j, i] = D_comp

            if compute_r:
                r_value = D_comp / (np.sqrt(p * (1 - p) * q * (1 - q)))
                r_matrix[i, j] = r_value
                r_matrix[j, i] = r_value

    np.fill_diagonal(ld_matrix, np.nan)

    if not compute_r:
        return ld_matrix
    else:
        return ld_matrix, r_matrix


def compute_LD_D_maladapt(
    gts: np.ndarray,
    ploidy: int = 1,
    compute_r: bool = True,
    maladapt_correction: bool = True,
) -> tuple[np.ndarray, np.ndarray] | np.ndarray:
    """
    Computes the linkage disequilibrium coefficient D and optionally the correlation coefficient r for pairs of SNPs,
    with an option to apply maladaptation correction.

    Parameters
    ----------
    gts : numpy.ndarray
        A 2D numpy array where each row represents a SNP and each column represents a genotype (haplotype) for an individual.
        The shape of the array is (n_snps, n_individuals), where `n_snps` is the number of SNPs and `n_individuals` is
        the number of individuals or haplotypes.

    ploidy : int, optional, default=1
        The ploidy level, indicating the number of chromosome sets. Typically, ploidy is 1 for haploids and 2 for diploids.

    compute_r : bool, optional, default=True
        If True, the function will compute the r coefficient, which measures the strength of correlation between the SNPs.
        If False, only the D coefficient is computed.

    maladapt_correction : bool, optional, default=True
        If True, the minor allele will be used for D in case the frequency of the major allele (1) is zero.

    Returns
    -------
    tuple of numpy.ndarray or numpy.ndarray
        - If `compute_r` is True, a tuple of two numpy arrays is returned:
          - `ld_matrix`: A symmetric matrix where entry (i, j) is the D value between SNP i and SNP j.
          - `r_matrix`: A symmetric matrix where entry (i, j) is the r value between SNP i and SNP j, or NaN if the r value
            cannot be computed.
        - If `compute_r` is False, only the `ld_matrix` is returned, which is a symmetric matrix with the D values.

    """
    num_snps, num_individuals = gts.shape

    # Initialize LD matrix with NaN
    ld_matrix = np.full((num_snps, num_snps), np.nan)
    r_matrix = np.zeros((num_snps, num_snps))

    for i in range(num_snps):
        for j in range(i + 1, num_snps):  # Compute only upper triangle (symmetry)
            geno1 = gts[i, :]
            geno2 = gts[j, :]

            P_11 = 0
            P_22 = 0
            P_12 = 0
            P_21 = 0

            for ind in range(num_individuals):
                g1 = geno1[ind]
                g2 = geno2[ind]

                if g1 == 0 and g2 == 0:
                    P_11 += 1
                elif g1 == ploidy and g2 == ploidy:
                    P_22 += 1
                elif g1 == ploidy and g2 < ploidy:
                    P_21 += 1
                elif g1 < ploidy and g2 == ploidy:
                    P_12 += 1

            P_11 /= num_individuals
            P_22 /= num_individuals
            P_12 /= num_individuals
            P_21 /= num_individuals

            P_1x = P_11 + P_12
            P_2x = P_22 + P_21
            P_x1 = P_11 + P_21
            P_x2 = P_12 + P_22

            D_comp = (P_11 * P_22) - (P_12 * P_21)

            ld_matrix[i, j] = D_comp
            ld_matrix[j, i] = D_comp

            if compute_r:
                denominator = np.sqrt(P_2x * (1 - P_2x) * P_x2 * (1 - P_x2))
                if denominator == 0:
                    if maladapt_correction:
                        denominator = np.sqrt(P_1x * (1 - P_1x) * P_x1 * (1 - P_x1))

                        r_value = D_comp / denominator
                    else:
                        r_value = np.nan
                else:
                    r_value = D_comp / denominator
                r_matrix[i, j] = r_value
                r_matrix[j, i] = r_value

    np.fill_diagonal(ld_matrix, np.nan)

    if compute_r:
        return ld_matrix, r_matrix
    else:
        return ld_matrix


# FILET features


def tajimas_d(gts: np.ndarray, ploidy: int = 1) -> float:
    """
    Compute Tajima's D statistic for a given genotype matrix.

    D = (θπ - θW) / sqrt(Var(θπ - θW))
    Formulae for parameters from Tajima 1989.

    Parameters
    ----------
    gts : np.ndarray
        A 2D NumPy array of shape (mutations x samples) where each entry represents the number of
        mutant alleles (0, 1, or 2) at each site for each individual.

    ploidy : int, optional, default=1
        The ploidy of the population. For example, ploidy = 2 for diploid organisms and ploidy = 1 for haploid organisms.

    Returns
    -------
    float
        Tajima's D statistic
    """
    theta_pi_val = theta_pi(gts, ploidy=ploidy)
    theta_W_val = theta_W(gts, ploidy=ploidy)
    S = num_segregating_sites(gts)

    num_snps, individuals = gts.shape
    a1 = sum(1 / i for i in range(1, individuals + 1))
    a2 = sum(1 / (i**2) for i in range(1, individuals + 1))

    b1 = (individuals + 1) / (3 * (individuals - 1))
    b2 = (2 * (individuals**2 + individuals + 3)) / (
        9 * individuals * (individuals - 1)
    )

    c1 = b1 - (1 / a1)
    c2 = b2 - ((individuals + 2) / (a1**2 * individuals)) + a2 / (a1**2)

    e1 = c1 / a1
    e2 = c2 / (a1**2 + a2)

    denominator = np.sqrt(e1 * S + e2 * S * (S - 1))

    if denominator == 0:
        return float("nan")

    result = (theta_pi_val - theta_W_val) / denominator

    return result


# additional features


def compute_D_plus(
    src_gts: np.ndarray,
    ref_gts: np.ndarray,
    tgt_gts: np.ndarray,
    out_gts: np.ndarray = None,
    ploidy: int = 1,
    compute_D_ancestral: bool = False,
) -> float:
    """
    Computes the D+-statistic for detecting admixture between populations (Lopez-Fang 2024).

    Parameters
    ----------
    src_gts : np.ndarray
        A 2D numpy array representing haplotypes of population 1 (source).
    ref_gts : np.ndarray
        A 2D numpy array representing haplotypes of population 2 (reference / sister group).
    tgt_gts : np.ndarray
        A 2D numpy array representing haplotypes of population 3 (target).
    out_gts : np.ndarray
        A 2D numpy array representing haplotypes of population 4 (outgroup).
        If not provided, it is assumed that the ancestral allel is always present in the outgroup, and thus the frequency of the derived allel (p4_freq) is 0.
    compute_D_ancestral : bool
        compute D_ancestral (which basically consists of the D+-terms w/o the standard D-terms, see Lopez-Fang 2024)

    Returns
    -------
    float
        The D+-statistic value (Lopez-Fang 2024).
    """

    # Compute allele frequencies using the provided calc_freq function
    src_freq = calc_freq(src_gts, ploidy=ploidy)
    ref_freq = calc_freq(ref_gts, ploidy=ploidy)
    tgt_freq = calc_freq(tgt_gts, ploidy=ploidy)
    if out_gts is not None:
        out_freq = calc_freq(out_gts, ploidy=ploidy)
    else:
        out_freq = 0

    # Compute ABBA and BABA site patterns
    abbavec = (1.0 - ref_freq) * tgt_freq * src_freq * (1.0 - out_freq)
    babavec = ref_freq * (1.0 - tgt_freq) * src_freq * (1.0 - out_freq)

    # Compute BAAA and ABAA site patterns
    baaavec = ref_freq * (1.0 - tgt_freq) * (1.0 - src_freq) * (1.0 - out_freq)
    abaavec = (1.0 - ref_freq) * tgt_freq * (1.0 - src_freq) * (1.0 - out_freq)

    # Summing up across loci
    abba = np.sum(abbavec)
    baba = np.sum(babavec)

    baaa = np.sum(baaavec)
    abaa = np.sum(abaavec)

    baaa_abaa_difference = baaa - abaa
    baaa_abaa_addition = baaa + abaa
    if compute_D_ancestral:
        return baaa_abaa_difference / baaa_abaa_addition

    abba_baba_difference = abba - baba
    abba_baba_addition = abba + baba

    D_plus = (abba_baba_difference + baaa_abaa_difference) / (
        abba_baba_addition + baaa_abaa_addition
    )
    return D_plus
