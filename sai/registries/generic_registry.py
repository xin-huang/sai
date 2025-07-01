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


from abc import ABC, abstractmethod
from typing import Any, Callable


class GenericRegistry(ABC):
    """
    Generic class for all registries.

    Provides a consistent interface for registering and retrieving
    named components such as missions, models, or features.
    """

    def __init__(self):
        self._registry = {}

    def register(self, name: str) -> Callable:
        """
        Decorator to register a class or callable under a given name.

        Parameters
        ----------
        name : str
            The key under which the component is registered.

        Returns
        -------
        Callable
            A decorator that registers the class/function.
        """

        def decorator(obj: Any) -> Any:
            self._register(name, obj)
            return obj

        return decorator

    def _register(self, name: str, obj: Any) -> None:
        if name in self._registry:
            raise ValueError(f"{name!r} is already registered.")
        self._registry[name] = obj

    def get(self, name: str) -> Any:
        """
        Retrieves a registered component by name.

        Parameters
        ----------
        name : str
            The key of the registered component.

        Returns
        -------
        Any
            The registered component.
        """
        if name not in self._registry:
            raise KeyError(f"No component registered under name '{name}'")
        return self._registry[name]

    def list_registered(self) -> list[str]:
        """
        Lists all registered component names.

        Returns
        -------
        list of str
            Names of all registered components.
        """
        return list(self._registry.keys())
