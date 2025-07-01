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


from pydantic import RootModel, field_validator
from typing import Dict, Union


class PloidyConfig(RootModel[Dict[str, Dict[str, int]]]):
    """
    Configuration for specifying per-population ploidy levels
    under categories like 'ref', 'tgt', 'src', and 'outgroup'.

    Ensures:
    - Only allowed keys are present
    - Each sub-dictionary maps to positive integers
    - Required keys ("ref", "tgt", "src") are present
    """

    @field_validator("root")
    def validate_ploidy_dict(
        cls, v: Dict[str, Dict[str, int]]
    ) -> Dict[str, Dict[str, int]]:
        allowed_keys = {"ref", "tgt", "src", "outgroup"}
        required_keys = {"ref", "tgt", "src"}

        extra_keys = set(v.keys()) - allowed_keys
        if extra_keys:
            raise ValueError(
                f"Unsupported ploidy keys: {extra_keys}. Allowed keys are {allowed_keys}."
            )

        missing_keys = required_keys - set(v.keys())
        if missing_keys:
            raise ValueError(f"Missing required ploidy keys: {missing_keys}.")

        for group, subdict in v.items():
            if not isinstance(subdict, dict):
                raise ValueError(
                    f"Value for '{group}' must be a dictionary of population -> ploidy."
                )
            for pop, ploidy in subdict.items():
                if not isinstance(ploidy, int) or ploidy <= 0:
                    raise ValueError(
                        f"Ploidy for '{group}:{pop}' must be a positive integer."
                    )

        return v

    def get_ploidy(self, group: str, population: str = None) -> Union[int, list[int]]:
        """
        Returns the ploidy for a given population under a given group.

        Parameters
        ----------
        group : str
            One of "ref", "tgt", "src", or "outgroup".
        population : str, optional
            The name of the population within the group. If None, return all ploidies as a list.

        Returns
        -------
        int or list[int]
            - If population is given: returns the ploidy for that population.
            - If population is None: returns a list of ploidies for all populations in the group.
        """
        if group not in self.root:
            raise KeyError(f"Group '{group}' not found in configuration.")

        if population is None:
            return list(self.root[group].values())

        if population not in self.root[group]:
            raise KeyError(
                f"Population '{population}' not found under group '{group}'."
            )

        return self.root[group][population]
