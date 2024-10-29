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


import os
import pandas as pd
from sai.utils.multiprocessing import mp_manager
from sai.utils.generators import WindowDataGenerator
from sai.utils.preprocessors import FeatureVectorsPreprocessor


def score(
    vcf_file: str,
    chr_name: str,
    ref_ind_file: str,
    tgt_ind_file: str,
    src_ind_file: str,
    win_len: int,
    win_step: int,
    num_src: int,
    anc_allele_file: str,
    ploidy: int,
    is_phased: bool,
    w: float,
    x: float,
    y: list[float],
    output_file: str,
    quantile: float,
    num_workers: int,
) -> None:
    """
    Processes and scores genomic data by generating windowed data and feature vectors.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file containing variant data.
    chr_name : str
        The chromosome name to be analyzed from the VCF file.
    ref_ind_file : str
        Path to the file containing reference population identifiers.
    tgt_ind_file : str
        Path to the file containing target population identifiers.
    src_ind_file : str
        Path to the file containing source population identifiers.
    win_len : int
        Length of each genomic window in base pairs.
    win_step : int
        Step size in base pairs between consecutive windows.
    num_src : int
        Number of source populations to include in each windowed analysis.
    anc_allele_file : str
        Path to the file containing ancestral allele information.
    ploidy : int
        The ploidy level of the genome.
    is_phased : bool
        Indicates if genotype data is phased.
    w : float
        Frequency threshold for calculating feature vectors.
    x : float
        Another frequency threshold for calculating feature vectors.
    y : list[float]
        List of frequency thresholds used for various calculations in feature vector processing.
    output_file : str
        File path to save the output of processed feature vectors.
    quantile : float
        Quantile threshold for feature vector processing.
    num_workers : int
        Number of parallel processes for multiprocessing.
    """
    generator = WindowDataGenerator(
        vcf_file=vcf_file,
        chr_name=chr_name,
        ref_ind_file=ref_ind_file,
        tgt_ind_file=tgt_ind_file,
        src_ind_file=src_ind_file,
        win_len=win_len,
        win_step=win_step,
        num_src=num_src,
        anc_allele_file=anc_allele_file,
        ploidy=ploidy,
        is_phased=is_phased,
    )

    preprocessor = FeatureVectorsPreprocessor(
        w=w,
        x=x,
        y=y,
        output_file=output_file,
        quantile=quantile,
    )

    header = f"Chrom\tStart\tEnd\tRef\tTgt\tSrc\tU\tQ{int(quantile*100)}\n"

    directory = os.path.dirname(output_file)
    if directory:
        os.makedirs(directory, exist_ok=True)
    with open(output_file, "w") as f:
        f.write(header)

    mp_manager(
        data_processor=preprocessor,
        data_generator=generator,
        nprocess=num_workers,
    )


def outlier(
    score_file: str, output_dir: str, output_prefix: str, quantile: float
) -> None:
    """
    Outputs rows exceeding the specified quantile for the second-to-last column (assumed to be 'U')
    and the last column (assumed to start with 'Q'), sorted by Start and then End columns.

    Parameters
    ----------
    score_file : str
        Path to the input file, in CSV format.
    output_dir : str
        Directory to store the output files.
    output_prefix : str
        Prefix for the output filenames.
    quantile : float
        Quantile threshold to filter rows.
    """
    # Read the input data file
    data = pd.read_csv(score_file, sep="\t")

    # Identify the U and Q columns as the second-to-last and last columns, respectively
    u_column = data.columns[-2]
    q_column = data.columns[-1]

    # Calculate quantile thresholds for the U and Q columns
    u_threshold = data[u_column].quantile(quantile)
    q_threshold = data[q_column].quantile(quantile)

    # Filter rows where values exceed the quantile thresholds
    u_outliers = data[data[u_column] > u_threshold]
    q_outliers = data[data[q_column] > q_threshold]

    # Sort the filtered data by 'Start' and then 'End' columns
    u_outliers_sorted = u_outliers.sort_values(by=["Start", "End"])
    q_outliers_sorted = q_outliers.sort_values(by=["Start", "End"])

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Define output file paths
    u_outliers_file = os.path.join(
        output_dir, f"{output_prefix}_{u_column}_outliers.tsv"
    )
    q_outliers_file = os.path.join(
        output_dir, f"{output_prefix}_{q_column}_outliers.tsv"
    )

    # Save the sorted filtered data to output files
    u_outliers_sorted.to_csv(u_outliers_file, index=False, sep="\t")
    q_outliers_sorted.to_csv(q_outliers_file, index=False, sep="\t")
