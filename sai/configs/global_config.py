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


from pydantic import BaseModel
from pydantic import model_validator
from typing import Optional
from sai.configs.stat_config import StatConfig
from sai.configs.ploidy_config import PloidyConfig
from sai.configs.pop_config import PopConfig


class GlobalConfig(BaseModel):
    statistics: StatConfig
    ploidies: PloidyConfig
    populations: PopConfig

    @model_validator(mode="after")
    def validate_population_in_ploidies(self) -> "GlobalConfig":
        """
        Cross-validates that every population in statistics also appears
        in the corresponding group in ploidies.
        """
        stat_data = (
            self.statistics.root
        )  # Dict[str, Dict[str, Dict[str, Union[float, tuple]]]]
        ploidy_data = self.ploidies.root  # Dict[str, Dict[str, int]]

        for stat_name, params in stat_data.items():
            for group in ("ref", "tgt", "src"):
                pop_dict = params.get(group, {})
                for pop in pop_dict:
                    if pop not in ploidy_data.get(group, {}):
                        raise ValueError(
                            f"Population '{pop}' used in statistics[{stat_name}][{group}] "
                            f"is not defined in ploidies[{group}]"
                        )
        return self

    @model_validator(mode="after")
    def validate_population_in_populations(self) -> "GlobalConfig":
        """
        Cross-validates that every population in statistics also appears
        in the corresponding group in sample files.
        """
        from sai.utils import parse_ind_file

        stat_data = self.statistics.root  # Dict[stat_name][group][pop] = ...
        population_paths = self.populations.root  # Dict[group] = path

        categories_per_group = {
            group: set(parse_ind_file(path).keys())
            for group, path in population_paths.items()
        }

        for stat_name, params in stat_data.items():
            for group in ("ref", "tgt", "src"):
                pop_dict = params.get(group, {})
                expected_categories = categories_per_group.get(group, set())

                for pop in pop_dict:
                    if pop not in expected_categories:
                        raise ValueError(
                            f"Population '{pop}' used in statistics[{stat_name}][{group}] "
                            f"is not found in the population file for group '{group}'."
                        )
        return self
