#!/usr/bin/env python

__author__ = "Daniel W. Davies"
__copyright__ = "Copyright Adam J. Jackson, Daniel W. Davies (2013)"
__version__ = "1.1"
__maintainer__ = "Daniel W. Davies"
__email__ = "d.w.davies@bath.ac.uk"
__date__ = "Jan 3 2017"

from setuptools import setup
import os

module_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    setup(
        name='SMACT',
        version='1.1.1',
        description='Semiconducting Materials by Analogy and Chemical Theory',
        long_description=open(os.path.join(module_dir, 'README.md')).read(),
        url='https://github.com/WMD-group/SMACT',
        author='Daniel W. Davies',
        author_email='d.w.davies@bath.ac.uk',
        license='GNU General Public License (GPL) v3',
        packages=['smact','smact.tests'],
        package_data={'smact': ['data/*.txt','data/*.csv']},
        zip_safe=False,
        install_requires=['scipy','numpy','spglib'],
        classifiers=['Programming Language :: Python :: 2.7',
                     'Development Status :: 5 - Production/Stable',
                     'Intended Audience :: Science/Research',
                     'Operating System :: OS Independent',
                     'Topic :: Scientific/Engineering']
    )
