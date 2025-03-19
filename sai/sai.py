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


import os
import pandas as pd
import matplotlib.pyplot as plt
from natsort import natsorted
from sai.utils.generators import ChunkGenerator
from sai.utils.preprocessors import ChunkPreprocessor


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
    w: float,
    x: float,
    y: list[float],
    output_file: str,
    stat_type: str,
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
    w : float
        Frequency threshold for calculating feature vectors.
    x : float
        Another frequency threshold for calculating feature vectors.
    y : list[float]
        List of frequency thresholds used for various calculations in feature vector processing.
    output_file : str
        File path to save the output of processed feature vectors.
    stat_type: str
        Specifies the type of statistic to compute.
    num_workers : int
        Number of parallel processes for multiprocessing.
    """
    generator = ChunkGenerator(
        vcf_file=vcf_file,
        chr_name=chr_name,
        window_size=win_len,
        step_size=win_step,
        num_chunks=num_workers * 8,
    )

    preprocessor = ChunkPreprocessor(
        vcf_file=vcf_file,
        ref_ind_file=ref_ind_file,
        tgt_ind_file=tgt_ind_file,
        src_ind_file=src_ind_file,
        win_len=win_len,
        win_step=win_step,
        w=w,
        x=x,
        y=y,
        output_file=output_file,
        stat_type=stat_type,
        anc_allele_file=anc_allele_file,
        num_src=num_src,
    )

    header = f"Chrom\tStart\tEnd\tRef\tTgt\tSrc\tnum_SNP\t{stat_type}\tCandidate\n"

    directory = os.path.dirname(output_file)
    if directory:
        os.makedirs(directory, exist_ok=True)
    with open(output_file, "w") as f:
        f.write(header)

    items = []

    for params in generator.get():
        items.extend(preprocessor.run(**params))

    preprocessor.process_items(items)


def outlier(score_file: str, output: str, quantile: float) -> None:
    """
    Outputs rows exceeding the specified quantile for the chosen column ('U' or 'Q'),
    sorted by Start and then End columns.

    Parameters
    ----------
    score_file : str
        Path to the input file, in CSV format.
    output : str
        Path to the output file.
    quantile : float
        Quantile threshold to filter rows.
    """
    # Read the input data file
    data = pd.read_csv(
        score_file, sep="\t", dtype={"Candidate Position": str}, index_col=False
    )

    column = data.columns[-2]

    # Convert column to numeric for computation
    data[column] = pd.to_numeric(data[column], errors="coerce")

    # Calculate quantile threshold for the chosen column
    threshold = data[column].quantile(quantile)

    # Filter rows where values exceed the quantile threshold
    outliers = data[data[column] > threshold]

    # Sort the filtered data by 'Chrom', 'Start', 'End' columns
    if not outliers.empty:
        outliers = outliers.reset_index(drop=True)
        outliers_sorted = outliers.iloc[
            natsorted(
                outliers.index,
                key=lambda i: (
                    outliers.loc[i, "Chrom"],
                    int(outliers.loc[i, "Start"]),
                    int(outliers.loc[i, "End"]),
                ),
            )
        ]
    else:
        outliers_sorted = outliers

    # Convert all columns to string before saving
    outliers_sorted = outliers_sorted.astype(str)

    # Save the sorted filtered data to the output file
    outliers_sorted.to_csv(output, index=False, sep="\t")


def plot(
    outlier_file: str,
    output: str,
    xlabel: str,
    ylabel: str,
    title: str,
    figsize_x: float = 6,
    figsize_y: float = 6,
    dpi: int = 300,
    alpha: float = 0.6,
) -> None:
    """
    Reads an outlier file and creates a scatter plot with U values on the Y-axis
    and Q values on the X-axis, then saves the plot to the specified output file.

    Parameters
    ----------
    outlier_file : str
        Path to the input file containing outlier data.
    output : str
        Path to save the output plot.
    xlabel : str
        Label for the X-axis.
    ylabel : str
        Label for the Y-axis.
    title : str
        Title of the plot.
    figsize_x : float, optional
        Width of the figure (default: 6).
    figsize_y : float, optional
        Height of the figure (default: 6).
    dpi : int, optional
        Resolution of the saved plot (default: 300).
    alpha : float, optional
        Transparency level of scatter points (default: 0.6).
    """
    # Read the input file
    data = pd.read_csv(outlier_file, sep="\t")

    # Identify the U and Q columns
    u_column = data.columns[-4]
    q_column = data.columns[-3]

    # Convert to numeric
    data[u_column] = pd.to_numeric(data[u_column], errors="coerce")
    data[q_column] = pd.to_numeric(data[q_column], errors="coerce")

    # Plot
    plt.figure(figsize=(figsize_x, figsize_y))
    plt.scatter(data[q_column], data[u_column], alpha=alpha)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(alpha=0.5, linestyle="--")

    # Save plot
    plt.savefig(output, dpi=dpi)
    plt.close()
