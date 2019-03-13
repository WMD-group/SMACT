[![DOI](https://zenodo.org/badge/14117740.svg)](https://zenodo.org/badge/latestdoi/14117740)
[![Documentation Status](https://readthedocs.org/projects/smact/badge/?version=latest)](http://smact.readthedocs.org/en/latest/?badge=latest)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://perso.crans.org/besson/LICENSE.html)
[![Build Status](https://travis-ci.org/WMD-group/SMACT.svg?branch=master)](https://travis-ci.org/WMD-group/SMACT)

SMACT
=====

**Semiconducting Materials from Analogy and Chemical Theory** (SMACT) is a collection of rapid screening tools that uses data about chemical elements.

![](SMACT.png)

*If you torture the data enough, nature will always confess* - Roland Coase (from 'How should economists choose?')

Statement of need
--------
The purpose of SMACT is to facilitate the high-throughput screening and design of functional materials. It follows a top-down approach where a set of element combinations is generated and then screened using rapid chemical filters. It can be used as part of a multi-technique workflow or to feed machine learning models for materials.

![](smact_simple.gif)

Contents
--------

* **smact** library containing:
  * **\_\_init\_\_.py** Contains the core `Element` and `Species` classes.
  *  **data_loader.py** Handles the loading of external data used to initialise the core `smact.Element` and `smact.Species` classes.
  *  **screening.py** Used for generating and applying filters to compositional search spaces.
  *  **properties.py** A collection of tools for estimating useful properties based on composition.
  * **lattice.py** Given the sites, multiplicities and possible oxidation states
    at those sites, this reads from the database and generates all possible
    stoichiometeries.
  * **builder.py** Builds some common lattice structures, given the chemical
    composition.
  * **lattice_parameters.py** Estimation of lattice parameters for various lattice types using covalent/ionic radii.
  * **distorter.py** A collection of functions for enumerating and then
    substituting on inequivalent sites of a sub-lattice.

Requirements
------------

The main language is Python 3 and basic requirements are Numpy and Scipy.
The [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase) (ASE),  [spglib](http://atztogo.github.io/spglib), and [pymatgen](www.pymatgen.org) are also required for many components.

Installation
------------
The latest stable release of SMACT can be installed via pip which will automatically setup other Python packages as required:

    pip install smact  

Alternatively, the very latest version can be installed using:

    pip install git+git://github.com/WMD-group/SMACT.git

For developer installation SMACT can be installed from a copy of the source
repository (https://github.com/wmd-group/smact); this will be preferred if using experimental code branches.

To clone the project from Github and make a local installation:

    git clone https://github.com/wmd-group/smact.git
    cd smact
    pip3 install --user -e .

With -e pip will create links to the source folder so that that changes
to the code will be immediately reflected on the PATH.

Usage
-----

SMACT's features are accessed through Python scripts, importing classes and functions as needed.
Some applications are available in [our examples folder](https://github.com/WMD-group/SMACT/tree/master/examples), as well as in the [SMACT workflow respository](https://github.com/WMD-group/SMACT_workflows).

License and attribution
-----------------------

Python code and original data tables are licensed under the GNU General Public License (GPL) v3.

The following files have their own licenses: **data/elements.txt** is from the [OpenBabel](http://openbabel.sourceforge.net) project and licensed under the GPL v2, which is included in the parent folder.

Development notes
-----------------

### Bugs, features and questions
Please use the [Issue Tracker](https://github.com/WMD-group/smact/issues) to report bugs or request features. While we hope that most questions can be answered by searching [the docs](https://smact.readthedocs.io/en/latest/), we welcome new questions on the issue tracker, especially if they helps us improve the docs!

### Code contributions

We are always looking for ways to make SMACT better and more useful to the wider community; contributions are very welcome. Please use the ["Fork and Pull"](https://guides.github.com/activities/forking/) workflow to make contributions and stick as closely as possible to the following:

- Code style should comply with [PEP8](http://www.python.org/dev/peps/pep-0008) where possible. [Google's house style](http://google-styleguide.googlecode.com/svn/trunk/pyguide.html)
is also helpful, including a good model for docstrings.
- Please use comments liberally when adding nontrivial features, and take the chance to clean up other people's code while looking at it.
- Add tests wherever possible, and use the test suite to check if you broke anything.
  Testing modules should be pass/fail and wrapped into **tests/test.py**.
  Run the tests using `python -m smact.tests.test -v`.
  (The final `-v` is optional and adds more detail to the output.)

References
----------

[D. W. Davies et al,
"Materials discovery by chemical analogy: role of oxidation states in structure prediction" *Faraday Discuss.* **211**, 553 (2018)](https://pubs.rsc.org/en/Content/ArticleLanding/2018/FD/C8FD00032H)

[D. W. Davies et al,
"Computer-aided design of metal chalcohalide semiconductors: from chemical composition to crystal structure" *Chem. Sci.* **9**, 1022 (2018)](http://www.cell.com/chem/abstract/S2451-9294(16)30155-3)

[D. W. Davies et al,
"Computational screening of all stoichiometric inorganic materials" *Chem* **1**, 617 (2016)](http://www.cell.com/chem/abstract/S2451-9294(16)30155-3)

[B. R. Pamplin, "A systematic method of deriving new semiconducting
compounds by structural analogy", *J. Phys. Chem. Solids*
**25**, 675 (1964)](http://www.sciencedirect.com/science/article/pii/0022369764901763)
