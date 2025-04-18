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


import os.path
from setuptools import setup, find_packages

# The directory containing this file
HERE = os.path.abspath(os.path.dirname(__file__))

# The text of the README file
with open(os.path.join(HERE, "README.md")) as fid:
    README = fid.read()

# This call to setup() does all the work
setup(
    name="sai-pg",
    python_requires="==3.9.19",
    version="1.0.1",
    description="A Python Package for Statistics for Adaptive Introgression",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/xin-huang/sai",
    author="Xin Huang",
    author_email="xinhuang.res@gmail.com",
    license="GPLv3",
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.9",
    ],
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "matplotlib==3.9.1",
        "natsort==8.4.0",
        "numpy==1.26.4",
        "pandas==2.2.1",
        "pysam==0.23.0",
        "scikit-allel==1.3.7",
        "scipy==1.12.0",
    ],
    entry_points={"console_scripts": ["sai=sai.__main__:main"]},
)
