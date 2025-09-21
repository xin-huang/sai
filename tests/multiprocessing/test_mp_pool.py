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


import pytest
from sai.multiprocessing.mp_pool import mp_pool, mp_worker


class DataPreprocessor:
    """Mock class simulating a data processor."""

    def run(self, x: int) -> int:
        """Mock processing: returns the square of x."""
        return x**2

    def process_items(self, results: list) -> None:
        """Stores the processed results for validation."""
        self.final_results = results


class DataGenerator:
    """Mock class simulating a data generator."""

    def __init__(self, data: list):
        self.data = data

    def get(self):
        """Yields data in dictionary format."""
        for x in self.data:
            yield {"x": x}


@pytest.mark.parametrize(
    "params, expected",
    [
        ((DataPreprocessor(), {"x": 2}), 4),
        ((DataPreprocessor(), {"x": 3}), 9),
        ((DataPreprocessor(), {"x": 4}), 16),
    ],
)
def test_mp_worker(params, expected):
    """Tests mp_worker to ensure correct processing."""
    assert mp_worker(params) == expected


def test_mp_pool():
    """Tests mp_pool to ensure parallel processing works correctly."""
    data_processor = DataPreprocessor()
    data_generator = DataGenerator([1, 2, 3, 4, 5])
    nprocess = 2

    # Run multiprocessing pool
    mp_pool(data_processor, data_generator, nprocess)

    # Validate results
    assert hasattr(data_processor, "final_results")  # Ensure results are stored
    assert sorted(data_processor.final_results) == [
        1,
        4,
        9,
        16,
        25,
    ]  # Check correctness
