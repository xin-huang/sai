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
import yaml
import pandas as pd
from pathlib import Path
from sai.generators import ChunkGenerator
from sai.preprocessors import ChunkPreprocessor
from sai.configs import GlobalConfig
from sai.utils.utils import natsorted_df


def score(
    vcf_file: str,
    chr_name: str,
    win_len: int,
    win_step: int,
    anc_allele_file: str,
    output_file: str,
    config: str,
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
    win_len : int
        Length of each genomic window in base pairs.
    win_step : int
        Step size in base pairs between consecutive windows.
    anc_allele_file : str
        Path to the file containing ancestral allele information.
    output_file : str
        File path to save the output of processed feature vectors.
    config: str
        Path to the YAML configuration file specifying the statistics and ploidies to compute.
    num_workers : int
        Number of parallel processes for multiprocessing.
    """
    try:
        with open(config, "r") as f:
            config_dict = yaml.safe_load(f)
    except FileNotFoundError:
        raise FileNotFoundError(f"Configuration file '{config}' not found.")
    except yaml.YAMLError as e:
        raise ValueError(f"Error parsing YAML configuration file '{config}': {e}")

    required_fields = ["statistics", "ploidies", "populations"]
    missing_fields = [field for field in required_fields if field not in config_dict]

    if missing_fields:
        raise ValueError(
            f"Missing required fields in configuration file '{config}': {', '.join(missing_fields)}"
        )

    global_config = GlobalConfig(**config_dict)

    stat_config = global_config.statistics
    ploidy_config = global_config.ploidies
    pop_config = global_config.populations

    generator = ChunkGenerator(
        vcf_file=vcf_file,
        chr_name=chr_name,
        window_size=win_len,
        step_size=win_step,
        # num_chunks=num_workers * 8,
        num_chunks=1,
    )

    preprocessor = ChunkPreprocessor(
        vcf_file=vcf_file,
        ref_ind_file=pop_config.get_population("ref"),
        tgt_ind_file=pop_config.get_population("tgt"),
        src_ind_file=pop_config.get_population("src"),
        out_ind_file=pop_config.get_population("outgroup"),
        win_len=win_len,
        win_step=win_step,
        output_file=output_file,
        ploidy_config=ploidy_config,
        stat_config=stat_config,
        anc_allele_file=anc_allele_file,
    )

    src_pops = list(ploidy_config.root["src"].keys())

    header_parts = ["Chrom", "Start", "End", "Ref", "Tgt", "Src", "N(Variants)"]

    for stat_name in stat_config.root.keys():
        if stat_name in ("U", "Q") or len(src_pops) <= 1:
            header_parts.append(stat_name)
        else:
            for sp in src_pops:
                header_parts.append(f"{stat_name}.{sp}")

    header = "\t".join(header_parts) + "\n"

    directory = os.path.dirname(output_file)
    if directory:
        os.makedirs(directory, exist_ok=True)
    with open(output_file, "w") as f:
        f.write(header)

    for key in ("U", "Q"):
        if key in stat_config.root:
            path = Path(output_file)
            log_file = path.with_suffix(f".{key}.log")
            with open(log_file, "w") as f:
                f.write(f"Chrom\tStart\tEnd\t{key}_SNP\n")

    items = []

    for params in generator.get():
        items.extend(preprocessor.run(**params))

    preprocessor.process_items(items)


def outlier(score_file: str, output_prefix: str, quantile: float) -> None:
    """
    Identifies outlier windows for each statistic column in a score file and
    write them to separate output files.

    This function reads a tab-delimited score file, determines which columns
    contain statistics (e.g., U, Q, D+, etc.), computes the specified quantile
    threshold for each statistic, and outputs rows exceeding that threshold.
    Results for each statistic are written to a separate TSV file, sorted by
    Chrom, Start, and End when available.

    Parameters
    ----------
    score_file : str
        Path to the input score file (tab-delimited).
    output_prefix : str
        Prefix for the output files. Each output file is named
        "{output_prefix}.{stat}.tsv".
    quantile : float
        Quantile threshold (between 0 and 1) used to define outliers.
    """
    df = pd.read_csv(score_file, sep="\t", na_values=["nan"], index_col=False)

    cols = list(df.columns)
    if "N(Variants)" in cols:
        start_idx = cols.index("N(Variants)") + 1
        metric_cols = cols[start_idx:]
    else:
        # fallback: exclude common non-metric columns, keep numeric ones
        non_metrics = {"Chrom", "Start", "End", "Ref", "Tgt", "Src"}
        candidate = [c for c in cols if c not in non_metrics]
        metric_cols = [
            c for c in candidate if pd.to_numeric(df[c], errors="coerce").notna().any()
        ]

    if not metric_cols:
        raise ValueError("No metric columns found.")

    for col in metric_cols:
        s_num = pd.to_numeric(df[col], errors="coerce").dropna()

        if s_num.empty:
            warnings.warn(
                f"Column '{col}' has no numeric values; writing empty result.",
                UserWarning,
            )
            out_sorted = pd.DataFrame(columns=df.columns)
        elif s_num.nunique() == 1:
            thr = s_num.iloc[0]
            warnings.warn(
                f"Column '{col}' has only one unique value ({thr}); writing empty result.",
                UserWarning,
            )
            out_sorted = pd.DataFrame(columns=df.columns)
        else:
            thr = s_num.quantile(quantile)
            col_num = pd.to_numeric(df[col], errors="coerce")
            if not col.startswith("U"):
                out = df[col_num >= thr]
            else:
                out = df[col_num > thr]

            if not out.empty:
                out = out.reset_index(drop=True)
                try:
                    out_sorted = natsorted_df(out)  # your existing natural sort
                except NameError:
                    keys = [k for k in ("Chrom", "Start", "End") if k in out.columns]
                    out_sorted = (
                        out.sort_values(by=keys, kind="mergesort") if keys else out
                    )
            else:
                out_sorted = out

        out_sorted.astype(str).to_csv(
            f"{output_prefix}.{col}.{quantile}.outliers.tsv", index=False, sep="\t"
        )
