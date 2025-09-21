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
from typing import Dict
from pydantic import RootModel, field_validator


REQUIRED_KEYS = {"ref", "tgt", "src"}
OPTIONAL_KEYS = {"outgroup"}
ALLOWED_KEYS = REQUIRED_KEYS | OPTIONAL_KEYS


class PopConfig(RootModel[Dict[str, str]]):
    """
    Configuration for population sample name files.

    Required:
        - ref: Path to file containing reference population sample names.
        - tgt: Path to file containing target population sample names.
        - src: Path to file containing source population sample names.

    Optional:
        - outgroup: Path to file containing outgroup sample names.
    """

    @field_validator("root")
    def validate_population_keys_and_paths(cls, v: Dict[str, str]) -> Dict[str, str]:
        keys = set(v.keys())
        missing = REQUIRED_KEYS - keys
        invalid = keys - ALLOWED_KEYS
        if missing:
            raise ValueError(f"Missing required population keys: {missing}")
        if invalid:
            raise ValueError(f"Unsupported population keys: {invalid}")
        for name, path in v.items():
            if not os.path.isfile(path):
                raise ValueError(f"{name} file does not exist: {path}")
        return v

    def get_population(self, group: str) -> str:
        """
        Retrieves the file path for a given population group.

        Parameters
        ----------
        group : str
            The population group name (e.g., 'ref', 'tgt', 'src', or 'outgroup').

        Returns
        -------
        str
            The file path corresponding to the group.

        Raises
        ------
        ValueError
            If the requested group is not present in the configuration.
        """
        if group not in self.root:
            if group == "outgroup":
                return None
            else:
                raise ValueError(f"Population group '{group}' not found in config.")
        return self.root[group]
