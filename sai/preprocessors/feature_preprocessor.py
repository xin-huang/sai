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
from pathlib import Path
from typing import Any
from sai.preprocessors import DataPreprocessor
from sai.registries.stat_registry import STAT_REGISTRY
from sai.configs import PloidyConfig, StatConfig


class FeaturePreprocessor(DataPreprocessor):
    """
    A preprocessor subclass for generating feature vectors from genomic data.

    This class extends DataPreprocessor to include additional functionality for creating
    feature vectors based on genomic variants, reference and target individual genotypes,
    and window-based genomic statistics.
    """

    def __init__(
        self,
        output_file: str,
        stat_config: StatConfig,
        anc_allele_available: bool = False,
    ):
        """
        Initializes FeatureVectorsPreprocessor with specific frequency thresholds
        and output file for storing generated feature vectors.

        Parameters
        ----------
        output_file : str
            Path to the output file to save processed feature vectors.
        stat_config: StatConfig,
            Specifies the configuration of statistics to compute.
        anc_allele_available: bool, optional
            If True, ancestral allele information is available.
            If False, ancestral allele information is unavailable.
            Default is False.
        """
        self.output_file = output_file
        self.anc_allele_available = anc_allele_available
        self.stat_config = stat_config

    def run(
        self,
        chr_name: str,
        ref_pop: str,
        tgt_pop: str,
        src_pop_list: list[str],
        start: int,
        end: int,
        pos: np.ndarray,
        ref_gts: np.ndarray,
        tgt_gts: np.ndarray,
        src_gts_list: list[np.ndarray],
        ploidy_config: PloidyConfig,
    ) -> list[dict[str, Any]]:
        """
        Generates feature vectors for a specified genomic window.

        Parameters
        ----------
        chr_name : str
            Chromosome name.
        ref_pop : str
            Reference population name.
        tgt_pop : str
            Target population name.
        src_pop_list : list[str]
            List of source population names.
        start : int
            Start position of the genomic window.
        end : int
            End position of the genomic window.
        pos : np.ndarray
            A 1D numpy array where each element represents the genomic position.
        ref_gts : np.ndarray
            Genotype data for the reference population.
        tgt_gts : np.ndarray
            Genotype data for the target population.
        src_gts_list : list[np.ndarray]
            List of genotype arrays for each source population.
        ploidy_config: PloidyConfig
            Configuration specifying ploidy levels for each population involved in the analysis.

        Returns
        -------
        list[dict[str, Any]]
            A list containing a dictionary of calculated feature vectors for the genomic window.
        """
        items = {
            "chr_name": chr_name,
            "start": start,
            "end": end,
            "ref_pop": ref_pop,
            "tgt_pop": tgt_pop,
            "src_pop_list": src_pop_list,
            "nsnps": len(pos),
            "cdd_pos": {},
        }

        if (
            (ref_gts is None)
            or (tgt_gts is None)
            or (src_gts_list is None)
            or (ploidy_config is None)
        ):
            for stat_name in self.stat_config.root.keys():
                if len(src_pop_list) > 1:
                    items[stat_name] = [np.nan for _ in range(len(src_pop_list))]
                else:
                    items[stat_name] = np.nan
                if stat_name in ["U", "Q"]:
                    items[stat_name] = np.nan
                    items["cdd_pos"][stat_name] = np.array([])
        else:
            for stat_name in self.stat_config.root.keys():
                stat_cls = STAT_REGISTRY.get(stat_name)
                stat = stat_cls(
                    ref_gts=ref_gts,
                    tgt_gts=tgt_gts,
                    src_gts_list=src_gts_list,
                    ref_ploidy=ploidy_config.get_ploidy("ref", ref_pop),
                    tgt_ploidy=ploidy_config.get_ploidy("tgt", tgt_pop),
                    src_ploidy_list=ploidy_config.get_ploidy("src"),
                )
                if stat_name == "U":
                    results = stat.compute(
                        pos=pos,
                        w=self.stat_config.get_parameters(stat_name)["ref"][ref_pop],
                        x=self.stat_config.get_parameters(stat_name)["tgt"][tgt_pop],
                        y_list=list(
                            self.stat_config.get_parameters(stat_name)["src"].values()
                        ),
                        anc_allele_available=self.anc_allele_available,
                    )
                    items["cdd_pos"][stat_name] = results["cdd_pos"]
                elif stat_name == "Q":
                    results = stat.compute(
                        pos=pos,
                        w=self.stat_config.get_parameters(stat_name)["ref"][ref_pop],
                        quantile=self.stat_config.get_parameters(stat_name)["tgt"][
                            tgt_pop
                        ],
                        y_list=list(
                            self.stat_config.get_parameters(stat_name)["src"].values()
                        ),
                        anc_allele_available=self.anc_allele_available,
                    )
                    items["cdd_pos"][stat_name] = results["cdd_pos"]
                else:
                    results = stat.compute()
                items[stat_name] = results["value"]

        return [items]

    def process_items(self, items: list[dict[str, Any]]) -> None:
        """
        Processes and writes a single dictionary of feature vectors to the output file.

        Parameters
        ----------
        items : dict[str, Any]
            A dictionary containing feature vectors for a genomic window.
        """
        with open(
            self.output_file, "a"
        ) as f:  # Open in append mode for continuous writing
            lines = []
            for item in items:
                src_pop = item["src_pop_list"]
                src_pop_str = ",".join(src_pop)

                stats_parts = []
                for stat_name in self.stat_config.root.keys():
                    val = item.get(stat_name)

                    if isinstance(val, list) and len(val) == len(src_pop):
                        # expand one value per src, in the same order as src_pop
                        stats_parts.extend("" if v is None else str(v) for v in val)
                    else:
                        # fallback: single value
                        if isinstance(val, list):
                            v = val[0] if len(val) > 0 else ""
                        else:
                            v = val
                        stats_parts.append("" if v is None else str(v))

                stats = "\t".join(stats_parts)

                line = (
                    f"{item['chr_name']}\t{item['start']}\t{item['end']}\t"
                    f"{item['ref_pop']}\t{item['tgt_pop']}\t{src_pop_str}\t"
                    f"{item['nsnps']}\t{stats}\n"
                )
                lines.append(line)

            f.writelines(lines)

        for key in ("U", "Q"):
            if key in self.stat_config.root:
                path = Path(self.output_file)
                log_file = path.with_suffix(f".{key}.log")
                with open(log_file, "a") as f:
                    lines = []
                    for item in items:
                        cdd = (
                            "NA"
                            if item["cdd_pos"][key].size == 0
                            else ",".join(
                                f"{item['chr_name']}:{pos}"
                                for pos in item["cdd_pos"][key]
                            )
                        )
                        line = f"{item['chr_name']}\t{item['start']}\t{item['end']}\t{cdd}\n"
                        lines.append(line)
                    f.writelines(lines)
