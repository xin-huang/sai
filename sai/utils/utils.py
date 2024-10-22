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
import numpy as np
from dataclasses import dataclass
from typing import Optional, Union


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
    anc_allele_file: Optional[str] = None,
    filter_missing: bool = True,
) -> dict[str, Union[ChromosomeData, list[str]]]:
    """
    Read genotype data from a VCF file efficiently, organized by chromosome.

    Parameters
    ----------
    vcf : str
        The name of the VCF file containing genotype data.
    ind_samples : dict of str to list of str
        A dictionary where keys are categories (e.g., different sample groups), and values are lists of sample names.
    anc_allele_file : str or None
        The name of the BED file containing ancestral allele information, or None if not provided.
    filter_missing : bool
        Indicates whether to filter out variants that are missing across all samples (default is True).

    Returns
    -------
    dict
        A dictionary containing genotype data organized by chromosome.
        Keys are chromosome names, each mapped to a ChromosomeData instance.
        Additionally, a "samples" key maps to a list of all sample names in the VCF.
    """
    try:
        # Load all samples from the VCF file
        all_samples = [sample for samples in ind_samples.values() for sample in samples]
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
        )
    except Exception as e:
        raise ValueError(f"Failed to read VCF file {vcf}: {e}")

    # Convert genotype data to a more efficient GenotypeArray
    gt = allel.GenotypeArray(vcf_data.get("calldata/GT"))
    chr_names = np.unique(vcf_data.get("variants/CHROM", []))
    pos = vcf_data.get("variants/POS")
    ref = vcf_data.get("variants/REF")
    alt = vcf_data.get("variants/ALT")

    if gt is None or pos is None or ref is None or alt is None:
        raise ValueError("Invalid VCF file: Missing essential genotype data fields.")

    # Load ancestral allele data if provided
    anc_allele = read_anc_allele(anc_allele_file) if anc_allele_file else None

    # Prepare indices grouped by chromosome for faster lookup
    chrom_indices = {c: np.where(vcf_data["variants/CHROM"] == c)[0] for c in chr_names}

    # Organize data by chromosome using ChromosomeData
    data = {}
    sample_indices = [all_samples.index(s) for s in all_samples]

    for c, indices in chrom_indices.items():
        # Create ChromosomeData for the chromosome
        chrom_data = ChromosomeData(
            POS=pos[indices],
            REF=ref[indices],
            ALT=alt[indices],
            GT=gt.take(indices, axis=0).take(sample_indices, axis=1),
        )

        # Remove missing data if specified
        if filter_missing:
            missing_index = chrom_data.GT.count_missing(axis=1) == len(sample_indices)
            chrom_data = filter_geno_data(chrom_data, ~missing_index)

        # Check and incorporate ancestral alleles if the file is provided
        if anc_allele:
            chrom_data = check_anc_allele(chrom_data, anc_allele, c)

        data[c] = chrom_data

    data["samples"] = vcf_data.get("samples")

    return data


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
    ref_ind_file: Optional[str],
    tgt_ind_file: Optional[str],
    src_ind_file: Optional[str],
    anc_allele_file: Optional[str],
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
]:
    """
    Helper function for reading data from reference, target, and source populations.

    Parameters
    ----------
    vcf_file : str
        Name of the VCF file containing genotype data.
    ref_ind_file : str or None
        File with reference population sample information. None if not provided.
    tgt_ind_file : str or None
        File with target population sample information. None if not provided.
    src_ind_file : str or None
        File with source population sample information. None if not provided.
    anc_allele_file : str or None
        File with ancestral allele information. None if not provided.
    is_phased : bool, optional
        Whether to use phased genotypes (default is True).
    filter_ref : bool, optional
        Whether to filter fixed variants for reference data (default is True).
    filter_tgt : bool, optional
        Whether to filter fixed variants for target data (default is True).
    filter_src : bool, optional
        Whether to filter fixed variants for source data (default is False).
    filter_missing : bool, optional
        Whether to filter out missing data (default is True).

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
        geno_data = read_geno_data(
            vcf_file, all_samples, anc_allele_file, filter_missing
        )
    except Exception as e:
        raise ValueError(f"Failed to read VCF data: {e}")

    # Separate reference, target, and source data
    ref_data = extract_group_data(geno_data, ref_samples)
    tgt_data = extract_group_data(geno_data, tgt_samples)
    src_data = extract_group_data(geno_data, src_samples)

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

    return ref_data, ref_samples, tgt_data, tgt_samples, src_data, src_samples


def extract_group_data(
    geno_data: dict[str, ChromosomeData],
    sample_groups: Optional[dict[str, list[str]]] = None,
) -> Optional[dict[str, dict[str, ChromosomeData]]]:
    """
    Extract genotype data from geno_data based on the sample groups.

    Parameters
    ----------
    geno_data : dict of str to ChromosomeData
        Contains genotype data organized by chromosome, with each value being a ChromosomeData instance.
    sample_groups : dict of str to list of str, optional
        Contains sample group information, where each key is a group name and the value is a list of samples.
        If None, the function returns None.

    Returns
    -------
    extracted_data : dict or None
        Genotype data organized by sample group and chromosome, or None if no sample groups are provided.
        The structure is as follows:

        - Outer keys represent sample group names.
        - Inner keys represent chromosome names.
        - Values are ChromosomeData instances, filtered to include only the samples in the specified group.
    """
    # If no sample groups are provided, return None
    if sample_groups is None:
        return None

    sample_indices = {sample: idx for idx, sample in enumerate(geno_data["samples"])}

    extracted_data = {}

    # Iterate over each sample group
    for group, samples in sample_groups.items():
        if group not in extracted_data:
            extracted_data[group] = {}

        # Initialize group data for each group
        group_data = {}

        # Iterate over each chromosome in geno_data
        for chrom, chrom_data in geno_data.items():
            # Skip the "samples" key in geno_data
            if chrom == "samples":
                continue

            # Extract sample indices
            indices = [sample_indices[s] for s in samples if s in sample_indices]

            # Extract genotype data for the selected samples
            group_data[chrom] = ChromosomeData(
                GT=chrom_data.GT[:, indices, :],
                POS=chrom_data.POS,
                REF=chrom_data.REF,
                ALT=chrom_data.ALT,
            )

        # Store the extracted data for the group
        extracted_data[group] = group_data

    return extracted_data


def filter_fixed_variants(
    data: dict[str, dict[str, ChromosomeData]], samples: dict[str, list[str]]
) -> dict[str, dict[str, ChromosomeData]]:
    """
    Filter out fixed variants for each population in the given data.

    Parameters
    ----------
    data : dict of str to dict of str to ChromosomeData
        Genotype data organized by category and chromosome, where each chromosome
        is represented by a ChromosomeData instance.
    samples : dict of str to list of str
        Sample information corresponding to each category, with each list containing
        sample names for a specific population category.

    Returns
    -------
    filtered_data : dict of str to dict of str to ChromosomeData
        Genotype data with fixed variants filtered out for each category and chromosome.
    """
    filtered_data = {}
    for cat, chrom_data in data.items():
        filtered_data[cat] = {}
        for c, geno in chrom_data.items():
            ref_fixed_variants = np.sum(geno.GT.is_hom_ref(), axis=1) == len(
                samples[cat]
            )
            alt_fixed_variants = np.sum(geno.GT.is_hom_alt(), axis=1) == len(
                samples[cat]
            )
            fixed_variants = np.logical_or(ref_fixed_variants, alt_fixed_variants)
            index = np.logical_not(fixed_variants)
            filtered_data[cat][c] = filter_geno_data(geno, index)

    return filtered_data


def reshape_genotypes(
    data: Optional[dict[str, dict[str, ChromosomeData]]], is_phased: bool
) -> None:
    """
    Reshape genotypes based on whether they are phased or unphased.

    Parameters
    ----------
    data : dict of str to dict of str to ChromosomeData or None
        Genotype data organized by sample group and chromosome. If None, the function does nothing.
    is_phased : bool
        If True, reshape phased genotypes. Otherwise, sum over ploidy for unphased genotypes.
    """
    if data is None:
        return

    for category in data:
        for c in data[category]:
            mut_num, ind_num, ploidy = data[category][c].GT.shape
            if is_phased:
                data[category][c].GT = np.reshape(
                    data[category][c].GT, (mut_num, ind_num * ploidy)
                )
            else:
                data[category][c].GT = np.sum(data[category][c].GT, axis=2)


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


def read_anc_allele(anc_allele_file: str) -> dict[str, dict[int, str]]:
    """
    Reads ancestral allele information from a BED file.

    Parameters
    ----------
    anc_allele_file : str
        Path to the BED file containing ancestral allele information.

    Returns
    -------
    dict of {str: dict of {int: str}}
        Chromosome-level dictionary mapping genomic positions to ancestral alleles.

    Raises
    ------
    FileNotFoundError
        If the ancestral allele file is not found.
    ValueError
        If no ancestral allele information is found in the file.
    """
    anc_allele = {}
    try:
        with open(anc_allele_file, "r") as f:
            for line in f:
                e = line.rstrip().split()
                chrom, pos, allele = e[0], int(e[2]), e[3]
                anc_allele.setdefault(chrom, {})[pos] = allele
    except FileNotFoundError:
        raise FileNotFoundError(f"File {anc_allele_file} not found.")

    if not anc_allele:
        raise ValueError("No ancestral allele is found! Please check your data.")

    return anc_allele


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
    chr_name: str,
    polymorphism_size: int,
    step_size: int = None,
    window_based: bool = True,
    random_polymorphisms: bool = False,
    seed: int = None,
) -> list[tuple]:
    """
    Creates sliding windows along the genome.

    Parameters
    ----------
    pos : np.ndarray
        Positions for the variants.
    chr_name : str
        Name of the chromosome.
    polymorphism_size : int
        Length of sliding windows or number of random positions.
    step_size : int, optional
        Step size of sliding windows. Default: None.
    window_based : bool, optional
        Whether to create sliding windows containing the start and end positions (True)
        or positions of each polymorphism within the window (False). Default: True.
    random_polymorphisms : bool, optional
        Whether to randomly select polymorphism positions (only used if window_based is False). Default: False.
    seed : int, optional
        Seed for the random number generator (only used if random_polymorphisms is True). Default: None.

    Returns
    -------
    list of tuple
        List of sliding windows along the genome if window_based is True,
        or list of position arrays if window_based is False. Each entry is either
        a tuple of (chr_name, start_position, end_position) for window_based, or
        (chr_name, numpy.ndarray of positions) for non-window_based.

    Raises
    ------
    ValueError
        If `step_size` or `polymorphism_size` are non-positive, or if the `pos` array is empty,
        or if no windows could be created with the given parameters.
    """
    if (step_size is not None and step_size <= 0) or polymorphism_size <= 0:
        raise ValueError(
            "`step_size` and `polymorphism_size` must be positive integers."
        )
    if len(pos) == 0:
        raise ValueError("`pos` array must not be empty.")

    window_positions = []

    if window_based:
        win_start = max(
            0, (pos[0] + step_size) // step_size * step_size - polymorphism_size
        )
        last_pos = pos[-1]

        while last_pos > win_start:
            win_end = win_start + polymorphism_size
            window_positions.append((chr_name, [win_start, win_end]))
            win_start += step_size
    else:
        if random_polymorphisms:
            if seed is not None:
                np.random.seed(seed)
            if len(pos) < polymorphism_size:
                raise ValueError(
                    "No windows could be created with the given number of polymorphisms."
                )
            polymorphism_indexes = np.random.choice(
                len(pos), size=polymorphism_size, replace=False
            )
            polymorphism_indexes = np.sort(polymorphism_indexes).tolist()
            window_positions.append((chr_name, polymorphism_indexes))
        else:
            i = 0
            while i + polymorphism_size <= len(pos):
                window_positions.append(
                    (chr_name, list(range(i, i + polymorphism_size)))
                )
                i += step_size

            if len(window_positions) == 0:
                raise ValueError(
                    "No windows could be created with the given number of polymorphisms and step size."
                )

            if window_positions[-1][1][1] != len(pos) - 1:
                window_positions.append(
                    (chr_name, list(range(len(pos) - polymorphism_size, len(pos))))
                )

    return window_positions
