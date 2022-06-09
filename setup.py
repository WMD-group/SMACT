#!/usr/bin/env python

__author__ = "Daniel W. Davies"
__copyright__ = "Copyright Daniel W. Davies, Adam J. Jackson, Keith T. Butler (2019)"
__version__ = "2.3.2"
__maintainer__ = "Anthony O. Onwuli"
__email__ = "anthony.onwuli16@imperial.ac.uk"
__date__ = "Nov 4 2021"

from setuptools import setup, Extension
import os
import unittest

module_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    setup(
        name="SMACT",
        version="2.3.2",
        description="Semiconducting Materials by Analogy and Chemical Theory",
        long_description=open(os.path.join(module_dir, "README.md")).read(),
        long_description_content_type="text/markdown",
        url="https://github.com/WMD-group/SMACT",
        author="Daniel W. Davies",
        author_email="d.davies16@imperial.ac.uk",
        maintainer="Anthony O. Onwuli",
        maintainer_email="anthony.onwuli16@imperial.ac.uk",
        license="MIT",
        packages=["smact", "smact.tests", "smact.structure_prediction"],
        package_data={
            "smact": [
                "data/*.txt",
                "data/*.csv",
                "data/*.data",
                "data/*.xlsx",
                "data/*.json",
            ]
        },
        zip_safe=False,
        test_suite="smact.tests.test",
        install_requires=[
            "scipy",
            "numpy",
            "spglib",
            "pymatgen",
            "ase",
            "pandas",
            "pathos",
        ],
        classifiers=[
            "Programming Language :: Python",
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Science/Research",
            "Operating System :: OS Independent",
            "Topic :: Scientific/Engineering",
        ],
    )
