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


import allel
import numpy as np
import pandas as pd
from typing import Optional, Union
from natsort import natsorted
from sai.utils.genomic_dataclasses import ChromosomeData


def parse_ind_file(filename: str) -> dict[str, list[str]]:
    """
    Read sample information from a file and organize it by categories.

    Parameters
    ----------
    filename : str
        The name of the file containing sample information.

    Returns
    -------
    samples : dict of str to list of str
        A dictionary where the keys represent categories, and the values are lists of samples
        associated with those categories.

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    ValueError
        If no samples are found in the file.
    """
    try:
        samples = {}

        with open(filename, "r") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) != 2:
                    continue

                category, sample = parts
                if category not in samples:
                    samples[category] = []
                samples[category].append(sample)

        if not samples:
            raise ValueError(f"No samples found in {filename}. Please check your data.")

    except FileNotFoundError:
        raise FileNotFoundError(
            f"File '{filename}' not found. Please check the file path."
        )

    return samples


def read_geno_data(
    vcf: str,
    ind_samples: dict[str, list[str]],
    chr_name: str,
    start: int = None,
    end: int = None,
    anc_allele_file: Optional[str] = None,
    filter_missing: bool = True,
) -> tuple[ChromosomeData, list[str], int]:
    """
    Read genotype data from a VCF file efficiently for a specified chromosome.

    Parameters
    ----------
    vcf : str
        The name of the VCF file containing genotype data.
    ind_samples : dict of str to list of str
        A dictionary where keys are categories (e.g., different sample groups), and values are lists of sample names.
    chr_name : str
        The name of the chromosome to read.
    start: int, optional
        The starting position (1-based, inclusive) on the chromosome. Default: None.
    end: int, optional
        The ending position (1-based, inclusive) on the chromosome. Default: None.
    anc_allele_file : str, optional
        The name of the BED file containing ancestral allele information, or None if not provided.
    filter_missing : bool, optional
        Indicates whether to filter out variants that are missing across all samples. Default: True.

    Returns
    -------
    chrom_data: ChromosomeData
        A ChromosomeData instance for the specified chromosome in the VCF.
    samples: list
        A list of samples in the data.
    ploidy: int
        Ploidy level of the organism.
    """
    try:
        # Load all samples from the VCF file
        all_samples = [sample for samples in ind_samples.values() for sample in samples]

        # Use region parameter to restrict to the specified chromosome
        if (start is None) and (end is None):
            region = f"{chr_name}"
        else:
            region = f"{chr_name}:{start}-{end}"
        vcf_data = allel.read_vcf(
            vcf,
            fields=[
                "calldata/GT",
                "variants/CHROM",
                "variants/POS",
                "variants/REF",
                "variants/ALT",
                "samples",
            ],
            alt_number=1,
            samples=all_samples,
            region=region,  # Specify the chromosome region
            tabix=None,
        )
    except Exception as e:
        raise ValueError(f"Failed to read VCF file {vcf} from {region}: {e}") from e

    # Convert genotype data to a more efficient GenotypeArray
    if vcf_data is None:
        return None, all_samples, None

    gt = allel.GenotypeArray(vcf_data.get("calldata/GT"))
    pos = vcf_data.get("variants/POS")
    ref = vcf_data.get("variants/REF")
    alt = vcf_data.get("variants/ALT")
    ploidy = gt.shape[2]

    if gt is None or pos is None or ref is None or alt is None:
        raise ValueError("Invalid VCF file: Missing essential genotype data fields.")

    # Load ancestral allele data if provided
    if anc_allele_file:
        anc_alleles = read_anc_allele(
            anc_allele_file=anc_allele_file,
            chr_name=chr_name,
            start=start,
            end=end,
        )
    else:
        anc_alleles = None

    sample_indices = [all_samples.index(s) for s in all_samples]

    chrom_data = ChromosomeData(
        POS=pos, REF=ref, ALT=alt, GT=gt.take(sample_indices, axis=1)
    )

    # Remove missing data if specified
    if filter_missing:
        non_missing_index = chrom_data.GT.count_missing(axis=1) == 0
        num_missing = len(non_missing_index) - np.sum(non_missing_index)
        if num_missing != 0:
            print(
                f"Found {num_missing} variants with missing genotypes, removing them ..."
            )
        chrom_data = filter_geno_data(chrom_data, non_missing_index)

    # Check and incorporate ancestral alleles if the file is provided
    if anc_alleles:
        chrom_data = check_anc_allele(chrom_data, anc_alleles, chr_name)

    return chrom_data, vcf_data.get("samples"), ploidy


def filter_geno_data(
    data: ChromosomeData, index: Union[np.ndarray, list[bool]]
) -> ChromosomeData:
    """
    Filter the genotype data based on the provided index.

    Parameters
    ----------
    data : ChromosomeData
        An instance of ChromosomeData containing genotype data, where each attribute corresponds to an array (e.g., POS, REF, ALT, GT).
    index : np.ndarray or list of bool
        A boolean or integer array indicating which rows to keep.

    Returns
    -------
    ChromosomeData
        A new ChromosomeData instance with filtered data, containing only the rows specified by the index.
    """
    return ChromosomeData(
        POS=data.POS[index],
        REF=data.REF[index],
        ALT=data.ALT[index],
        GT=data.GT.compress(index, axis=0),
    )


def read_data(
    vcf_file: str,
    chr_name: str,
    ref_ind_file: Optional[str],
    tgt_ind_file: Optional[str],
    src_ind_file: Optional[str],
    anc_allele_file: Optional[str],
    start: int = None,
    end: int = None,
    is_phased: bool = True,
    filter_ref: bool = True,
    filter_tgt: bool = True,
    filter_src: bool = False,
    filter_missing: bool = True,
) -> tuple[
    Optional[dict[str, dict[str, ChromosomeData]]],
    Optional[dict[str, list[str]]],
    Optional[dict[str, dict[str, ChromosomeData]]],
    Optional[dict[str, list[str]]],
    Optional[dict[str, dict[str, ChromosomeData]]],
    Optional[dict[str, list[str]]],
    Optional[int],
]:
    """
    Helper function for reading data from reference, target, and source populations.

    Parameters
    ----------
    vcf_file : str
        Name of the VCF file containing genotype data.
    chr_name: str
        Name of the chromosome to read.
    ref_ind_file : str or None
        File with reference population sample information. None if not provided.
    tgt_ind_file : str or None
        File with target population sample information. None if not provided.
    src_ind_file : str or None
        File with source population sample information. None if not provided.
    anc_allele_file : str or None
        File with ancestral allele information. None if not provided.
    start: int, optional
        The starting position (1-based, inclusive) on the chromosome. Default: None.
    end: int, optional
        The ending position (1-based, inclusive) on the chromosome. Default: None.
    is_phased : bool, optional
        Whether to use phased genotypes. Default: True.
    filter_ref : bool, optional
        Whether to filter fixed variants for reference data. Default: True.
    filter_tgt : bool, optional
        Whether to filter fixed variants for target data. Default: True.
    filter_src : bool, optional
        Whether to filter fixed variants for source data. Default: False.
    filter_missing : bool, optional
        Whether to filter out missing data. Default: True.

    Returns
    -------
    ref_data : dict or None
        Genotype data from reference populations, organized by category and chromosome.
    ref_samples : dict or None
        Sample information from reference populations.
    tgt_data : dict or None
        Genotype data from target populations, organized by category and chromosome.
    tgt_samples : dict or None
        Sample information from target populations.
    src_data : dict or None
        Genotype data from source populations, organized by category and chromosome.
    src_samples : dict or None
        Sample information from source populations.
    ploidy: int or None
        Ploidy level of the organism.

    Notes
    -----
    The `ref_data`, `tgt_data`, and `src_data` are organized as nested dictionaries where:

        - The outermost keys represent different populations or sample categories.
        - The second-level keys represent different chromosomes.
        - The innermost value is a ChromosomeData instance containing:
            - "POS": numpy array of variant positions.
            - "REF": numpy array of reference alleles.
            - "ALT": numpy array of alternative alleles.
            - "GT": allel.GenotypeArray containing genotype data.

    This organization allows easy access and manipulation of genotype data by category and chromosome,
    enabling flexible processing across different populations and chromosomal regions.
    """
    ref_data = ref_samples = tgt_data = tgt_samples = src_data = src_samples = None

    # Parse sample information
    if ref_ind_file:
        ref_samples = parse_ind_file(ref_ind_file)

    if tgt_ind_file:
        tgt_samples = parse_ind_file(tgt_ind_file)

    if src_ind_file:
        src_samples = parse_ind_file(src_ind_file)

    # Combine all samples for a single VCF read
    all_samples = {}
    if ref_samples:
        all_samples.update(ref_samples)
    if tgt_samples:
        all_samples.update(tgt_samples)
    if src_samples:
        all_samples.update(src_samples)

    try:
        # Read VCF data
        geno_data, all_samples, ploidy = read_geno_data(
            vcf=vcf_file,
            ind_samples=all_samples,
            chr_name=chr_name,
            start=start,
            end=end,
            anc_allele_file=anc_allele_file,
            filter_missing=filter_missing,
        )
    except Exception as e:
        raise ValueError(f"Failed to read VCF data: {e}")

    if geno_data is None:
        return None, ref_samples, None, tgt_samples, None, src_samples, None

    # Separate reference, target, and source data
    ref_data = extract_group_data(geno_data, all_samples, ref_samples)
    tgt_data = extract_group_data(geno_data, all_samples, tgt_samples)
    src_data = extract_group_data(geno_data, all_samples, src_samples)

    # Apply fixed variant filtering conditionally
    if filter_ref and ref_data and ref_samples:
        ref_data = filter_fixed_variants(ref_data, ref_samples)
    if filter_tgt and tgt_data and tgt_samples:
        tgt_data = filter_fixed_variants(tgt_data, tgt_samples)
    if filter_src and src_data and src_samples:
        src_data = filter_fixed_variants(src_data, src_samples)

    # Adjust genotypes based on phased/unphased requirement
    reshape_genotypes(ref_data, is_phased)
    reshape_genotypes(tgt_data, is_phased)
    reshape_genotypes(src_data, is_phased)

    return ref_data, ref_samples, tgt_data, tgt_samples, src_data, src_samples, ploidy


def extract_group_data(
    geno_data: dict[str, ChromosomeData],
    all_samples: list[str],
    sample_groups: Optional[dict[str, list[str]]] = None,
) -> Optional[dict[str, ChromosomeData]]:
    """
    Extract genotype data from geno_data based on the sample groups.

    Parameters
    ----------
    geno_data : dict of str to ChromosomeData
        Contains genotype data, where each value is a ChromosomeData instance.
    all_samples: list
        A list of all sample names in the dataset.
    sample_groups : dict of str to list of str, optional
        Contains sample group information, where each key is a group name and the value is a list of samples.
        If None, the function returns None.

    Returns
    -------
    extracted_data : dict or None
        Genotype data organized by sample group, or None if no sample groups are provided.
        The structure is as follows:

        - Keys represent sample group names.
        - Values are ChromosomeData instances, filtered to include only the samples in the specified group.
    """
    if sample_groups is None:
        return None

    sample_indices = {sample: idx for idx, sample in enumerate(all_samples)}

    extracted_data = {}

    for group, samples in sample_groups.items():
        indices = [sample_indices[s] for s in samples if s in sample_indices]

        # Extract ChromosomeData for the selected samples in each group
        extracted_data[group] = ChromosomeData(
            GT=geno_data.GT[:, indices, :],
            POS=geno_data.POS,
            REF=geno_data.REF,
            ALT=geno_data.ALT,
        )

    return extracted_data


def filter_fixed_variants(
    data: dict[str, ChromosomeData], samples: dict[str, list[str]]
) -> dict[str, ChromosomeData]:
    """
    Filter out fixed variants for each population in the given data.

    Parameters
    ----------
    data : dict of str to ChromosomeData
        Genotype data organized by category, where each category is represented by a ChromosomeData instance.
    samples : dict of str to list of str
        Sample information corresponding to each category, with each list containing
        sample names for a specific population category.

    Returns
    -------
    filtered_data : dict of str to ChromosomeData
        Genotype data with fixed variants filtered out for each category.
    """
    filtered_data = {}
    for cat, geno in data.items():
        ref_fixed_variants = np.sum(geno.GT.is_hom_ref(), axis=1) == len(samples[cat])
        alt_fixed_variants = np.sum(geno.GT.is_hom_alt(), axis=1) == len(samples[cat])
        fixed_variants = np.logical_or(ref_fixed_variants, alt_fixed_variants)
        index = np.logical_not(fixed_variants)
        filtered_data[cat] = filter_geno_data(geno, index)

    return filtered_data


def reshape_genotypes(
    data: Optional[dict[str, ChromosomeData]], is_phased: bool
) -> None:
    """
    Reshape genotypes based on whether they are phased or unphased.

    Parameters
    ----------
    data : dict of str to ChromosomeData or None
        Genotype data organized by sample group. If None, the function does nothing.
    is_phased : bool
        If True, reshape phased genotypes. Otherwise, sum over ploidy for unphased genotypes.
    """
    if data is None:
        return

    for category, chrom_data in data.items():
        mut_num, ind_num, ploidy = chrom_data.GT.shape
        if is_phased:
            chrom_data.GT = np.reshape(chrom_data.GT, (mut_num, ind_num * ploidy))
        else:
            chrom_data.GT = np.sum(chrom_data.GT, axis=2)


def get_ref_alt_allele(
    ref: list[str], alt: list[str], pos: list[int]
) -> tuple[dict[int, str], dict[int, str]]:
    """
    Indexes REF and ALT alleles with genomic positions.

    Parameters
    ----------
    ref : list of str
        REF alleles.
    alt : list of str
        ALT alleles.
    pos : list of int
        Genomic positions.

    Returns
    -------
    Dictionaries mapping genomic positions to REF and ALT alleles, respectively.
    """
    return {p: r for p, r in zip(pos, ref)}, {p: a for p, a in zip(pos, alt)}


def read_anc_allele(
    anc_allele_file: str, chr_name: str, start: int = None, end: int = None
) -> dict[str, dict[int, str]]:
    """
    Reads ancestral allele information from a BED file for a specified chromosome,
    optionally within a specified position range.

    Parameters
    ----------
    anc_allele_file : str
        Path to the BED file containing ancestral allele information.
    chr_name : str
        Name of the chromosome to read.
    start : int, optional
        Start position (1-based, inclusive) of the region to filter. If None, no lower bound.
    end : int, optional
        End position (1-based, inclusive) of the region to filter. If None, no upper bound.

    Returns
    -------
    dict of {str: dict of {int: str}}
        Chromosome-level dictionary mapping genomic positions to ancestral alleles.

    Raises
    ------
    FileNotFoundError
        If the ancestral allele file is not found.
    ValueError
        If no ancestral allele information is found for the specified chromosome (and region if specified).
    """
    anc_alleles = {}
    try:
        with open(anc_allele_file, "r") as f:
            for line in f:
                e = line.rstrip().split()
                chrom, pos, allele = e[0], int(e[2]), e[3]
                if chrom != chr_name:
                    continue
                if (start is not None and pos < start) or (
                    end is not None and pos > end
                ):
                    continue
                anc_alleles.setdefault(chrom, {})[pos] = allele
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"File {anc_allele_file} not found.") from exc

    if not anc_alleles:
        if start is not None or end is not None:
            raise ValueError(
                f"No ancestral allele is found for chromosome {chr_name} in the region {start}-{end}."
            )
        else:
            raise ValueError(f"No ancestral allele is found for chromosome {chr_name}.")

    return anc_alleles


def check_anc_allele(
    data: dict[str, ChromosomeData], anc_allele: dict[str, dict[int, str]], c: str
) -> dict[str, ChromosomeData]:
    """
    Checks whether the REF or ALT allele is the ancestral allele and updates genotype data.

    Parameters
    ----------
    data : dict
        Genotype data for checking ancestral allele information.
    anc_allele : dict of {str: dict of {int: str}}
        Dictionary with ancestral allele information.
    c : str
        Chromosome name.

    Returns
    -------
    dict
        Genotype data with updated alleles after ancestral allele checking.
    """
    ref_allele, alt_allele = get_ref_alt_allele(data.REF, data.ALT, data.POS)

    # Determine variants to remove or flip
    intersect_snps = np.intersect1d(
        list(ref_allele.keys()), list(anc_allele.get(c, {}).keys())
    )
    removed_snps, flipped_snps = [], []

    for v in intersect_snps:
        if anc_allele[c][v] not in {ref_allele[v], alt_allele[v]}:
            removed_snps.append(v)
        elif anc_allele[c][v] == alt_allele[v]:
            flipped_snps.append(v)

    # Filter data by intersecting SNPs and remove any that should be removed
    intersect_filter = np.in1d(data.POS, intersect_snps)
    data = filter_geno_data(data, intersect_filter)

    if removed_snps:
        remain_filter = np.logical_not(np.in1d(data.POS, removed_snps))
        data = filter_geno_data(data, remain_filter)

    # Flip alleles in SNPs where ALT allele is ancestral
    flip_snps(data, flipped_snps)

    return data


def flip_snps(data: dict[str, ChromosomeData], flipped_snps: list[int]) -> None:
    """
    Flips the genotypes for SNPs where the ALT allele is the ancestral allele.

    Parameters
    ----------
    data : dict
        Genotype data.
    flipped_snps : list of int
        List of positions where the ALT allele is ancestral.
    """
    # Create a boolean mask for positions that need to be flipped
    is_flipped = np.isin(data.POS, flipped_snps)

    # Flip all genotypes at once where the mask is True
    data.GT[is_flipped] = allel.GenotypeArray(abs(data.GT[is_flipped] - 1))


def split_genome(
    pos: np.ndarray,
    window_size: int,
    step_size: int,
    start: int = None,
) -> list[tuple]:
    """
    Creates sliding windows along the genome based on variant positions.

    Parameters
    ----------
    pos : np.ndarray
        Array of positions for the variants.
    window_size : int
        Length of each sliding window.
    step_size : int
        Step size of the sliding windows.
    start: int, optional
        Minimum starting coordinate for the first window. The first window will start
        no smaller than this value. Default is None.

    Returns
    -------
    list of tuple
        List of sliding windows, where each entry is a tuple (start_position, end_position)
        representing the start and end positions of each window.

    Raises
    ------
    ValueError
        - If `step_size` or `window_size` are non-positive
        - If `step_size` is greater than `window_size`
        - If the `pos` array is empty
    """
    # Validate inputs
    if step_size <= 0 or window_size <= 0:
        raise ValueError("`step_size` and `window_size` must be positive integers.")
    if step_size > window_size:
        raise ValueError("`step_size` cannot be greater than `window_size`.")
    if len(pos) == 0:
        raise ValueError("`pos` array must not be empty.")

    window_positions = []
    win_start = (pos[0] + step_size) // step_size * step_size - window_size + 1
    if start is None:
        start = 1
    win_start = max(win_start, start)

    # Create windows based on step size and window size
    while win_start <= pos[-1]:
        win_end = win_start + window_size - 1
        window_positions.append((win_start, win_end))
        win_start += step_size

    return window_positions


def natsorted_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sorts a DataFrame naturally by "Chrom", "Start", and "End" columns.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to be sorted.

    Returns
    -------
    pd.DataFrame
        The naturally sorted DataFrame.

    Raises
    ------
    ValueError
        If the required columns "Chrom", "Start", or "End" are missing.
    """
    required_columns = {"Chrom", "Start", "End"}

    if missing_columns := required_columns - set(df.columns):
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")

    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)

    sorted_indices = natsorted(
        df.index, key=lambda i: (df.at[i, "Chrom"], df.at[i, "Start"], df.at[i, "End"])
    )

    return df.loc[sorted_indices].reset_index(drop=True)
