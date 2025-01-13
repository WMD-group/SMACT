#!/usr/bin/env python

"""Installation for SMACT."""

from __future__ import annotations

__author__ = "The SMACT Developers"
__author_email__ = "a.walsh@imperial.ac.uk"
__copyright__ = "Copyright Daniel W. Davies, Adam J. Jackson, Keith T. Butler (2019)"
__version__ = "3.0.2"
__maintainer__ = "Anthony O. Onwuli"
__maintainer_email__ = "anthony.onwuli16@imperial.ac.uk"
__date__ = "January 13 2025"


import os

from setuptools import setup

module_dir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(module_dir, "README.md")) as f:
    __long_description__ = f.read()

if __name__ == "__main__":
    setup(
        name="SMACT",
        version=__version__,
        description="Semiconducting Materials by Analogy and Chemical Theory",
        long_description=__long_description__,
        long_description_content_type="text/markdown",
        url="https://github.com/WMD-group/SMACT",
        author=__author__,
        author_email=__author_email__,
        maintainer=__maintainer__,
        maintainer_email=__maintainer_email__,
        license="MIT",
        packages=[
            "smact",
            "smact.utils",
            "smact.tests",
            "smact.structure_prediction",
            "smact.dopant_prediction",
            "smact.utils.crystal_space",
        ],
        package_data={
            "smact": [
                "data/*.txt",
                "data/*.csv",
                "data/*.data",
                "data/*.xlsx",
                "data/*.json",
                "data/species_rep/*.json",
            ]
        },
        zip_safe=False,
        test_suite="smact.tests.test",
        install_requires=[
            "scipy",
            "numpy<3",
            "spglib",
            "pymatgen>=2024.2.20",
            "ase",
            "pandas",
            "pathos",
            "typing-extensions",
        ],
        classifiers=[
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.10",
            "Programming Language :: Python :: 3.11",
            "Programming Language :: Python :: 3.12",
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Science/Research",
            "Operating System :: OS Independent",
            "License :: OSI Approved :: MIT License",
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Chemistry",
        ],
        python_requires=">=3.10",
    )
