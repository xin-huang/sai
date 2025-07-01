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
import warnings
import numpy as np
import pandas as pd
from natsort import natsorted
from typing import Optional, Union
from sai.utils.genomic_dataclasses import ChromosomeData
from sai.configs import PloidyConfig


def parse_ind_file(filename: str) -> dict[str, list[str]]:
    """
    Read sample information from a file and organize it by categories.

    Parameters
    ----------
    filename : str
        The name of the file containing sample information.

    Returns
    -------
    samples : dict[str, list[str]]
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
    ploidy: int = 2,
    start: int = None,
    end: int = None,
    anc_allele_file: Optional[str] = None,
    filter_missing: bool = True,
) -> dict[str, ChromosomeData]:
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
    ploidy : int, optional
        Ploidy level of the genome.
    start : int, optional
        The starting position (1-based, inclusive) on the chromosome. Default: None.
    end : int, optional
        The ending position (1-based, inclusive) on the chromosome. Default: None.
    anc_allele_file : str, optional
        The name of the BED file containing ancestral allele information, or None if not provided.
    filter_missing : bool, optional
        Indicates whether to filter out variants that are missing across all samples. Default: True.

    Returns
    -------
    A dictionary mapping each population name to its ChromosomeData.
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
            numbers={"GT": ploidy},
            region=region,  # Specify the chromosome region
            tabix=None,
        )
    except Exception as e:
        raise ValueError(f"Failed to read VCF file {vcf} from {region}: {e}") from e

    if vcf_data is None:
        return None

    gt = allel.GenotypeArray(vcf_data.get("calldata/GT"))
    pos = vcf_data.get("variants/POS")
    ref = vcf_data.get("variants/REF")
    alt = vcf_data.get("variants/ALT")
    sample_names = list(vcf_data.get("samples"))

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

    chrom_data_dict = {}

    for pop, pop_samples in ind_samples.items():
        indices = [sample_names.index(s) for s in pop_samples]
        pop_gt = gt.take(indices, axis=1)

        chrom_data = ChromosomeData(
            POS=pos.copy(), REF=ref.copy(), ALT=alt.copy(), GT=pop_gt
        )

        missing_mask = chrom_data.GT.count_missing(axis=1) != 0

        if filter_missing:
            if np.any(missing_mask):
                chrom_data = filter_geno_data(chrom_data, ~missing_mask)
        else:
            if np.any(missing_mask):
                raise ValueError(
                    "Missing data is found. Please remove variants with missing data or enable filtering."
                )

        if anc_alleles:
            chrom_data = check_anc_allele(chrom_data, anc_alleles, chr_name)

        chrom_data_dict[pop] = chrom_data

    return chrom_data_dict


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
    ploidy_config: PloidyConfig,
    ref_ind_file: Optional[str],
    tgt_ind_file: Optional[str],
    src_ind_file: Optional[str],
    out_ind_file: Optional[str],
    anc_allele_file: Optional[str],
    start: int = None,
    end: int = None,
    is_phased: bool = True,
    filter_ref: bool = True,
    filter_tgt: bool = True,
    filter_src: bool = False,
    filter_out: bool = False,
    filter_missing: bool = True,
) -> dict[
    str,
    tuple[
        Optional[dict[str, dict[str, ChromosomeData]]],
        Optional[dict[str, list[str]]],
        Optional[dict[str, dict[str, ChromosomeData]]],
        Optional[dict[str, list[str]]],
        Optional[dict[str, dict[str, ChromosomeData]]],
        Optional[dict[str, list[str]]],
    ],
]:
    """
    Helper function for reading data from reference, target, and source populations.

    Parameters
    ----------
    vcf_file : str
        Name of the VCF file containing genotype data.
    chr_name : str
        Name of the chromosome to read.
    ploidy_config : PloidyConfig
        Configuration specifying ploidy levels for each population involved in the analysis.
    ref_ind_file : str or None
        File with reference population sample information. None if not provided.
    tgt_ind_file : str or None
        File with target population sample information. None if not provided.
    src_ind_file : str or None
        File with source population sample information. None if not provided.
    out_ind_file : str or None
        File with outgroup population sample information. None if not provided.
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
    filter_out : bool, optional
        Whether to filter fixed variants for outgroup data. Default: False.
    filter_missing : bool, optional
        Whether to filter out missing data. Default: True.

    Returns
    -------
    result : dict
        {
            "ref": (ref_data, ref_samples),
            "tgt": (tgt_data, tgt_samples),
            "src": (src_data, src_samples),
            "outgroup": (out_data, out_samples)  # optional
        }

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
    out_data : dict or None
        Genotype data from outgroup populations, organized by category and chromosome.
    out_samples : dict or None
        Sample information from outgroup populations.

    Notes
    -----
    The `ref_data`, `tgt_data`, `src_data`, `out_data` are organized as nested dictionaries where:

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
    group_params = [
        ("ref", ref_ind_file, filter_ref),
        ("tgt", tgt_ind_file, filter_tgt),
        ("src", src_ind_file, filter_src),
        ("outgroup", out_ind_file, filter_out),
    ]

    results = {}

    for group, ind_file, filter_flag in group_params:
        if ind_file is None:
            results[group] = (None, None)
            continue

        if group == "outgroup" and group not in ploidy_config.root:
            results[group] = (None, None)
            continue

        data, samples = _load_population_data(
            vcf_file=vcf_file,
            chr_name=chr_name,
            sample_file=ind_file,
            anc_allele_file=anc_allele_file,
            start=start,
            end=end,
            is_phased=is_phased,
            filter_flag=filter_flag,
            filter_missing=filter_missing,
            ploidy_config=ploidy_config,
            group=group,
        )
        results[group] = (data, samples)

    return results


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


def _load_population_data(
    vcf_file: str,
    chr_name: str,
    sample_file: Optional[str],
    anc_allele_file: Optional[str],
    start: Optional[int],
    end: Optional[int],
    is_phased: bool,
    filter_flag: bool,
    filter_missing: bool,
    ploidy_config: PloidyConfig,
    group: str,  # e.g., "ref", "tgt", "src"
) -> tuple[
    Optional[dict[str, dict[str, ChromosomeData]]], Optional[dict[str, list[str]]]
]:
    """
    Loads genotype data and sample information for a population group (e.g., reference) from a VCF file,
    handling multiple populations with potentially different ploidy.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file containing variant data.
    chr_name : str
        Chromosome name to extract from the VCF.
    sample_file : str or None
        Path to the file containing sample IDs grouped by population.
        If None, no data is loaded.
    anc_allele_file : str or None
        Path to the BED file with ancestral allele annotations.
    start : int or None
        Start position on the chromosome (1-based, inclusive). If None, starts at the beginning.
    end : int or None
        End position on the chromosome (1-based, inclusive). If None, reads to the end.
    is_phased : bool
        Whether the genotypes are phased.
    filter_flag : bool
        Whether to remove variants fixed in all samples of each population.
    filter_missing : bool
        Whether to filter out variants with missing genotypes across all samples.
    ploidy_config : PloidyConfig
        Configuration containing ploidy for all populations in all groups.
    group : str
        The group label (e.g., "ref", "tgt", "src") used to extract populations from ploidy_config.

    Returns
    -------
    data : dict[str, dict[str, ChromosomeData]] or None
        Dictionary mapping population -> chromosome -> ChromosomeData.
    samples : dict[str, list[str]] or None
        Dictionary mapping population -> list of sample IDs.
    """
    if sample_file is None:
        return None, None

    samples = parse_ind_file(sample_file)

    if group not in ploidy_config.root:
        raise ValueError(f"Ploidy configuration missing group '{group}'.")

    group_ploidies = ploidy_config.root[group]

    # Ensure all populations in ploidy_config[group] are in sample_file
    for population in group_ploidies:
        if population not in samples:
            raise ValueError(
                f"Population '{population}' in ploidy_config[{group}] not found in sample file: {sample_file}"
            )

    data: dict[str, ChromosomeData] = {}

    for population, sample_list in samples.items():
        if population not in group_ploidies:
            warnings.warn(
                f"Population '{population}' found in sample file but not in ploidy_config[{group}]; skipping.",
                RuntimeWarning,
            )
            continue

        ploidy = group_ploidies[population]

        try:
            geno_data = read_geno_data(
                vcf=vcf_file,
                ind_samples={population: sample_list},
                chr_name=chr_name,
                start=start,
                end=end,
                anc_allele_file=anc_allele_file,
                filter_missing=filter_missing,
                ploidy=ploidy,
            )
        except Exception as e:
            raise ValueError(
                f"Failed to read VCF data for {sample_file}, population '{population}': {e}"
            )

        if geno_data is None:
            continue

        if filter_flag:
            geno_data = filter_fixed_variants(geno_data, {population: sample_list})

        reshape_genotypes(geno_data, is_phased)

        data[population] = geno_data[
            population
        ]  # geno_data: dict[population -> ChromosomeData]

    if not data:
        return None, samples

    return data, samples
