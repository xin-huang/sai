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


from sai.multiprocessing import mp_manager
from sai.preprocessors import DataPreprocessor
from sai.generators import DataGenerator


class TmpDataPreprocessor(DataPreprocessor):
    def run(self, rep):
        return [rep]

    def process_items(self, items):
        print(items)


class TmpDataGenerator(DataGenerator):
    def __init__(self, start_rep=0, nrep=5):
        self.start_rep = start_rep
        self.nrep = nrep

    def get(self):
        for i in range(self.start_rep, self.start_rep + self.nrep):
            yield {"rep": i}

    def __len__(self):
        return self.nrep


class FailureDataPreprocessor(DataPreprocessor):
    def run(self, rep):
        raise Exception("Simulating failure by stopping.")

    def process_items(self, items):
        print(items)


def test_mp_manager(capfd):
    nprocess = 2
    nrep = 5

    data_processor = TmpDataPreprocessor()
    generator = TmpDataGenerator(nrep=nrep)

    mp_manager(
        data_processor=data_processor, data_generator=generator, nprocess=nprocess
    )

    # Define the expected set of outputs
    expected_set = {"[0, 1, 2, 3, 4]"}

    # Capture the actual output and convert it to a set of strings
    captured = capfd.readouterr()
    actual_set = {captured.out.strip()}

    # Compare the actual and expected sets
    assert actual_set == expected_set, "The output does not match the expected results."


def test_mp_manager_failure(capfd):
    nprocess = 2
    data_processor = FailureDataPreprocessor()
    generator = TmpDataGenerator(nrep=5)

    mp_manager(
        data_processor=data_processor, data_generator=generator, nprocess=nprocess
    )

    # Use capfd to capture stdout and stderr
    captured = capfd.readouterr()

    # Assertions to verify expected output and behavior
    assert "Simulating failure by stopping." in captured.err
    assert "did not complete successfully. Initiating shutdown." in captured.out
    assert "All workers are terminated." in captured.out
