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
import warnings
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from sai.utils.generators import ChunkGenerator
from sai.utils.preprocessors import ChunkPreprocessor
from sai.utils.utils import natsorted_df


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

    header = f"Chrom\tStart\tEnd\tRef\tTgt\tSrc\tN(SNP)\t{stat_type}\tCandidate\n"

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
        score_file,
        sep="\t",
        na_values=["nan"],
        dtype={"Candidate": str},
        index_col=False,
    )

    column = data.columns[-2]

    # Convert column to numeric for computation
    data[column] = pd.to_numeric(data[column], errors="coerce")

    # Calculate quantile threshold for the chosen column
    threshold = data[column].quantile(quantile)

    if data[column].nunique() == 1:
        warnings.warn(
            f"Column '{column}' contains only one unique value ({threshold}), making quantile filtering meaningless.",
            UserWarning,
        )
        outliers = pd.DataFrame(columns=data.columns)
    elif (threshold == 1) and (column.startswith("Q")):
        outliers = data[data[column] >= threshold]
    else:
        outliers = data[data[column] > threshold]

    # Sort the filtered data by 'Chrom', 'Start', 'End' columns
    if not outliers.empty:
        outliers = outliers.reset_index(drop=True)
        outliers_sorted = natsorted_df(outliers)
    else:
        outliers_sorted = outliers

    # Convert all columns to string before saving
    outliers_sorted = outliers_sorted.astype(str)

    # Save the sorted filtered data to the output file
    outliers_sorted.to_csv(output, index=False, sep="\t")


def plot(
    u_file: str,
    q_file: str,
    output: str,
    xlabel: str,
    ylabel: str,
    title: str,
    figsize_x: float = 6,
    figsize_y: float = 6,
    dpi: int = 300,
    alpha: float = 0.6,
    marker_size: float = 20,
    marker_color: str = "blue",
    marker_style: str = "o",
) -> None:
    """
    Reads two outlier files (U and Q), finds common candidate positions, and plots U vs. Q.

    Parameters
    ----------
    u_file : str
        Path to the input file containing U score/outlier data.
    q_file : str
        Path to the input file containing Q score/outlier data.
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
    marker_size : float, optional
        Size of the scatter plot markers (default: 20). See matplotlib.pyplot.scatter.
    marker_color : str, optional
        Color of the markers (default: "blue"). See matplotlib.pyplot.scatter.
    marker_style : str, optional
        Shape of the marker, e.g., 'o' for circle, '^' for triangle, 's' for square (default: "o").
        See matplotlib.pyplot.scatter for available marker styles.
    """

    # Read input files
    u_data = pd.read_csv(u_file, sep="\t")
    q_data = pd.read_csv(q_file, sep="\t")

    # Identify columns
    u_column = u_data.columns[-2]
    q_column = q_data.columns[-2]

    # Create genomic interval key as `chr:start-end`
    u_data["interval"] = (
        u_data["Chrom"].astype(str)
        + ":"
        + u_data["Start"].astype(str)
        + "-"
        + u_data["End"].astype(str)
    )
    q_data["interval"] = (
        q_data["Chrom"].astype(str)
        + ":"
        + q_data["Start"].astype(str)
        + "-"
        + q_data["End"].astype(str)
    )

    # Convert U and Q columns to numeric
    u_data[u_column] = pd.to_numeric(u_data[u_column], errors="coerce")
    q_data[q_column] = pd.to_numeric(q_data[q_column], errors="coerce")
    u_data = u_data.dropna(subset=[u_column])
    q_data = q_data.dropna(subset=[q_column])

    # Convert to dictionary: {interval: value}
    u_interval_dict = {row["interval"]: row[u_column] for _, row in u_data.iterrows()}
    q_interval_dict = {row["interval"]: row[q_column] for _, row in q_data.iterrows()}

    # Find intersection of intervals
    common_intervals = set(u_interval_dict.keys()) & set(q_interval_dict.keys())

    if not common_intervals:
        raise ValueError(
            "No common genomic intervals found between U and Q outlier files. The plot will be empty.",
        )

    # Create DataFrame for intersection and save to file
    intersection_df = pd.DataFrame(
        {
            "Chrom": [interval.split(":")[0] for interval in common_intervals],
            "Start": [
                int(interval.split(":")[1].split("-")[0])
                for interval in common_intervals
            ],
            "End": [
                int(interval.split(":")[1].split("-")[1])
                for interval in common_intervals
            ],
            u_column: [u_interval_dict[c] for c in common_intervals],
            q_column: [q_interval_dict[c] for c in common_intervals],
        }
    )

    # Save intersection data
    intersection_df_sorted = natsorted_df(intersection_df)
    intersection_output = os.path.splitext(output)[0] + ".intersection.tsv"
    intersection_df_sorted.to_csv(intersection_output, sep="\t", index=False)

    # Plot
    plt.figure(figsize=(figsize_x, figsize_y))
    plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.scatter(
        x=intersection_df[q_column],
        y=intersection_df[u_column],
        alpha=alpha,
        s=marker_size,
        c=marker_color,
        marker=marker_style,
    )
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(alpha=alpha, linestyle="--")

    # Save plot
    plt.savefig(output, dpi=dpi)
    plt.close()
