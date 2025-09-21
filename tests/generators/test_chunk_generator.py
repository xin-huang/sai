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
from sai.generators import ChunkGenerator


def test_chunk_generator():
    # Initialize
    generator = ChunkGenerator(
        vcf_file="tests/data/test.data.vcf",
        chr_name="21",
        step_size=5000,
        window_size=10000,
        num_chunks=2,
    )

    # Check that length is calculated properly (mocked to 3 records)
    assert len(generator) == 2  # num_workers

    # Check that chunks were split correctly
    expected_chunks = [(1, 30000), (25001, 55000)]
    assert generator.chunks == expected_chunks


def test_chunk_generator_chr_not_found():
    with pytest.raises(ValueError, match="Chromosome 1 not found in VCF."):
        ChunkGenerator(
            vcf_file="tests/data/test.data.vcf",
            chr_name="1",
            step_size=10000,
            window_size=10000,
            num_chunks=2,
        )
