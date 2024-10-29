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


import scipy
import yaml
import numpy as np
from typing import Any
from sai.stats.features import calc_rd
from sai.stats.features import calc_u
from sai.stats.features import calc_q
from sai.utils import parse_ind_file
from sai.utils.preprocessors import DataPreprocessor


class FeatureVectorsPreprocessor(DataPreprocessor):
    """
    A preprocessor subclass for generating feature vectors from genomic data.

    This class extends DataPreprocessor to include additional functionality for creating
    feature vectors based on genomic variants, reference and target individual genotypes,
    and window-based genomic statistics.

    """

    def __init__(self, ref_ind_file: str, tgt_ind_file: str, feature_config: str):
        """
        Initializes a new instance of FeatureVectorsPreprocessor with specific parameters.

        Parameters:
        -----------
        ref_ind_file : str
            Path to the file listing reference individual identifiers.
        tgt_ind_file : str
            Path to the file listing target individual identifiers.
        feature_config : str
            Path to the configuration file specifying the features to be computed.

        Raises
        ------
        FileNotFoundError
            If the feature configuration file is not found.
        ValueError
            If the feature configuration file is incorrectly formatted or does not contain any features.

        """
        try:
            with open(feature_config, "r") as f:
                features = yaml.safe_load(f)
            self.features = features.get("Features", {})
            if not self.features:
                raise ValueError("No features found in the configuration.")
        except FileNotFoundError:
            raise FileNotFoundError(
                f"Feature configuration file {feature_config} not found."
            )
        except yaml.YAMLError as exc:
            raise ValueError(f"Error parsing feature configuration: {exc}")

        ref_samples = parse_ind_file(ref_ind_file)
        tgt_samples = parse_ind_file(tgt_ind_file)
        self.samples = {
            "Ref": ref_samples,
            "Tgt": tgt_samples,
        }

    def run(
        self,
        chr_name: str,
        start: int,
        end: int,
        ploidy: int,
        is_phased: bool,
        ref_gts: np.ndarray,
        tgt_gts: np.ndarray,
        pos: np.ndarray,
    ) -> list[dict[str, Any]]:
        """
        Executes the feature vector generation process for a specified genomic window.

        Parameters
        ----------
        chr_name : str
            Name of the chromosome.
        start : int
            Start position of the genomic window.
        end : int
            End position of the genomic window.
        ploidy : int
            Ploidy of the samples, typically 2 for diploid organisms.
        is_phased : bool
            Indicates whether the genomic data is phased.
        ref_gts : np.ndarray
            Genotype array for the reference individuals.
        tgt_gts : np.ndarray
            Genotype array for the target individuals.
        pos : np.ndarray
            Array of variant positions within the genomic window.

        Returns
        -------
        list
            A list of dictionaries containing the formatted feature vectors for the genomic window.

        """
        variants_not_in_ref = np.sum(ref_gts, axis=1) == 0
        sub_ref_gts = ref_gts[variants_not_in_ref]
        sub_tgt_gts = tgt_gts[variants_not_in_ref]
        sub_pos = pos[variants_not_in_ref]

        items = dict()
        items["chr_name"] = chr_name
        items["start"] = start
        items["end"] = end
        if self.features.get("Total variant number", False):
            items["Total_var_num"] = cal_mut_num(ref_gts, tgt_gts, mut_type="total")

        if self.features.get("Private variant number", False):
            items["Private_var_num"] = cal_mut_num(
                sub_ref_gts, sub_tgt_gts, mut_type="private"
            )

        if self.features.get("Spectrum", False):
            spectra = cal_n_ton(tgt_gts, is_phased=is_phased, ploidy=ploidy)
            items["Spectrum"] = spectra

        if ("Ref distances" in self.features) and (self.features["Ref distances"]):
            ref_dists = cal_dist(ref_gts, tgt_gts)
            items.update(self._cal_dist_stats(ref_dists, "Ref"))

        if ("Tgt distances" in self.features) and (self.features["Tgt distances"]):
            tgt_dists = cal_dist(tgt_gts, tgt_gts)
            items.update(self._cal_dist_stats(tgt_dists, "Tgt"))

        if ("Sstar" in self.features) and (self.features["Sstar"]):
            sstar_scores, sstar_snp_nums, haplotypes = cal_sstar(
                sub_tgt_gts,
                sub_pos,
                method=self.features["Sstar"]["Genotype distance"],
                match_bonus=self.features["Sstar"]["Match bonus"],
                max_mismatch=self.features["Sstar"]["Max mismatch"],
                mismatch_penalty=self.features["Sstar"]["Mismatch penalty"],
            )
            items["Sstar"] = sstar_scores

        fmtted_res = self._fmt_res(
            res=items,
            ploidy=ploidy,
            is_phased=is_phased,
            samples=self.samples,
            features=self.features,
        )

        return fmtted_res

    def _cal_dist_stats(self, dists: np.ndarray, pop: str) -> dict:
        """
        Calculates statistical metrics for distance arrays.

        Parameters
        ----------
        dists : np.ndarray
            The distance array for which to calculate statistics.
        pop : str
            The population identifier ('ref' for reference, 'tgt' for target)
            used as a prefix in the results dictionary.

        Returns
        -------
        dict
            A dictionary of calculated statistical metrics, keyed by metric name.

        """
        stats = {}
        if self.features[f"{pop} distances"].get("All", False):
            stats[f"All_{pop}_dists"] = dists
        if self.features[f"{pop} distances"].get("Minimum", False):
            stats[f"Minimum_{pop}_dists"] = np.min(dists, axis=1)
        if self.features[f"{pop} distances"].get("Maximum", False):
            stats[f"Maximum_{pop}_dists"] = np.max(dists, axis=1)
        if self.features[f"{pop} distances"].get("Mean", False):
            stats[f"Mean_{pop}_dists"] = np.mean(dists, axis=1)
        if self.features[f"{pop} distances"].get("Median", False):
            stats[f"Median_{pop}_dists"] = np.median(dists, axis=1)
        if self.features[f"{pop} distances"].get("Variance", False):
            stats[f"Variance_{pop}_dists"] = np.var(dists, axis=1)
        if self.features[f"{pop} distances"].get("Skew", False):
            stats[f"Skew_{pop}_dists"] = scipy.stats.skew(dists, axis=1)
            stats[f"Skew_{pop}_dists"][np.isnan(stats[f"Skew_{pop}_dists"])] = 0
        if self.features[f"{pop} distances"].get("Kurtosis", False):
            stats[f"Kurtosis_{pop}_dists"] = scipy.stats.kurtosis(dists, axis=1)
            stats[f"Kurtosis_{pop}_dists"][np.isnan(stats[f"Kurtosis_{pop}_dists"])] = 0

        return stats

    def _fmt_res(
        self,
        res: dict[str, Any],
        ploidy: int,
        is_phased: bool,
        samples: dict[list[str]],
        features: dict[str, Any],
    ) -> list[dict[str, Any]]:
        """
        Formats the result dictionaries into a pandas DataFrame with appropriate headers.

        Parameters
        ----------
        res : dict[str, Any]
            A dictionary representing the results for a genomic window.
        ploidy : int
            The ploidy of the samples being processed.
        is_phased : bool
            Indicates whether the genomic data is phased.
        samples : dict[list[str]]
            A dictionary of reference and target sample identifiers.
        features : dict[str, Any]
            A dictionary specifying which features were calculated and should be included in the output.

        Returns
        -------
        list
            A list of dictionaries containing the formatted results with one row per sample and one column per feature.

        """
        if is_phased:
            num_samples = len(samples["Tgt"]) * ploidy
        else:
            num_samples = len(samples["Tgt"])

        base_dict = {
            "Chromosome": res["chr_name"],
            "Start": res["start"],
            "End": res["end"],
        }
        sample_dicts = [base_dict.copy() for _ in range(num_samples)]
        for i, sample_dict in enumerate(sample_dicts):
            sample = (
                f'{samples["Tgt"][int(i/ploidy)]}_{i%ploidy+1}'
                if is_phased
                else samples["Tgt"][i]
            )
            sample_dict["Sample"] = sample

            if ("Sstar" in features) and (features["Sstar"]):
                sample_dict["Sstar"] = res["Sstar"][i]
            if features.get("Total variant number", False):
                sample_dict["Total_var_num"] = res["Total_var_num"][i]
            if features.get("Private variant number", False):
                sample_dict["Private_var_num"] = res["Private_var_num"][i]
            if features.get("Spectrum", False):
                for j in range(num_samples + 1):
                    sample_dict[f"{j}_ton"] = res["Spectrum"][i][j]

            for pop in ["Ref", "Tgt"]:
                if (f"{pop} distances" in features) and (features[f"{pop} distances"]):
                    sample_dict.update(
                        self._fmt_dist_res(
                            row=res,
                            idx=i,
                            is_phased=is_phased,
                            ploidy=ploidy,
                            samples=self.samples,
                            features=self.features,
                            pop=pop,
                        )
                    )

        return sample_dicts

    def _fmt_dist_res(
        self,
        row: dict[str, Any],
        idx: int,
        is_phased: bool,
        ploidy: int,
        samples: dict[list[str]],
        features: dict[str, Any],
        pop: str,
    ) -> dict[str, float]:
        """
        Formats distance-related results for a single sample and population (reference or target).

        Parameters
        ----------
        row : dict[str, Any]
            A dictionary containing results for a genomic window.
        idx : int
            The index of the sample in the target list.
        is_phased : bool
            Indicates whether the genomic data is phased.
        ploidy : int
            The ploidy of the samples being processed.
        samples : dict[list[str]]
            A dictionary of reference and target sample identifiers.
        features : dict[str, Any]
            A dictionary specifying which features were calculated and should be included in the output.
        pop : str
            Indicates the population type ('Ref' or 'Tgt') for which distances are being formatted.

        Returns
        -------
        dict[str, float]
            A dictionary with keys for each distance-related feature and values for the current sample.

        """
        items = dict()

        dist_features = features.get(f"{pop} distances", {})
        for feature in [
            "Minimum",
            "Maximum",
            "Mean",
            "Median",
            "Variance",
            "Skew",
            "Kurtosis",
        ]:
            if dist_features.get(feature, False):
                items[f"{feature}_{pop}_dist"] = row[f"{feature}_{pop}_dists"][idx]

        if dist_features.get("All", False):
            if is_phased:
                for j in range(len(row[f"All_{pop}_dists"][idx])):
                    sample = samples[pop][int(j / ploidy)]
                    items[f"{pop}_dist_{sample}_{j%ploidy+1}"] = row[
                        f"All_{pop}_dists"
                    ][idx][j]
            else:
                for j in range(len(row[f"All_{pop}_dists"][idx])):
                    sample = samples[pop][j]
                    items[f"{pop}_dist_{sample}"] = row[f"All_{pop}_dists"][idx][j]

        return items
