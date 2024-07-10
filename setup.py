#!/usr/bin/env python

__author__ = "Daniel W. Davies"
__author_email__ = "d.w.davies@imperial.ac.uk"
__copyright__ = (
    "Copyright Daniel W. Davies, Adam J. Jackson, Keith T. Butler (2019)"
)
__version__ = "2.6"
__maintainer__ = "Anthony O. Onwuli"
__maintaier_email__ = "anthony.onwuli16@imperial.ac.uk"
__date__ = "July 10 2024"

import os
import unittest

from setuptools import Extension, setup

module_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    setup(
        name="SMACT",
        version=__version__,
        description="Semiconducting Materials by Analogy and Chemical Theory",
        long_description=open(os.path.join(module_dir, "README.md")).read(),
        long_description_content_type="text/markdown",
        url="https://github.com/WMD-group/SMACT",
        author=__author__,
        author_email=__author_email__,
        maintainer=__maintainer__,
        maintainer_email=__maintaier_email__,
        license="MIT",
        packages=[
            "smact",
            "smact.tests",
            "smact.structure_prediction",
            "smact.dopant_prediction",
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
            "numpy<2",
            "spglib",
            "pymatgen>=2024.2.20",
            "ase",
            "pandas",
            "pathos",
            "typing-extensions",
        ],
        classifiers=[
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.9",
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
        python_requires=">=3.9",
    )
