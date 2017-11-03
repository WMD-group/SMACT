SMACT
=====

**Semiconducting Materials from Analogy and Chemical Theory** (SMACT) is a collection of fast screening tools from elemental data.

![](SMACT.png)

*If you torture the data enough, nature will always confess* - Roland Coase (from 'How should economists choose?')

[![Documentation Status](https://readthedocs.org/projects/smact/badge/?version=latest)](http://smact.readthedocs.org/en/latest/?badge=latest)

[![DOI](https://zenodo.org/badge/14117740.svg)](https://zenodo.org/badge/latestdoi/14117740)

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

The main language is Python 3 with Numpy, Scipy and Matplotlib.
The [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase) 
(ASE) is required for some components, as is [spglib](http://atztogo.github.io/spglib).

The [chemlab](http://chemlab.github.com/chemlab) project is not
currently used, but is considered "friendly"; we will try to avoid
namespace clashes and it may be used for some features in the future.
Needed modifications will be made in [ajjackson's
fork](https://github.com/ajjackson/chemlab), 
but are expected to make it upstream fairly rapidly.

Installation
------------
	pip install git+git://github.com/WMD-group/SMACT.git

On a unix-like system, simply add the directory containing this README file
to your PYTHONPATH. e.g. in ~/.bashrc

    export PYTHONPATH="/home/username/src/smact:$PYTHONPATH"

Usage
-----

SMACT's features are
accessed through Python scripts, importing classes and functions as needed.
Some practical applications using are available in [our examples folder](https://github.com/WMD-group/SMACT/tree/master/examples).

License and attribution
-----------------------

SMACT is produced by the Walsh Materials Design group. Python code
and original data tables are licensed under the GNU General Public
License (GPL) v3.

The following files have their own licenses: **data/elements.txt** is
from the [OpenBabel](http://openbabel.sourceforge.net) project and licensed under the GPL v2, which is
included in the parent folder.

References
----------

[D. W. Davies et al, 
"Computational Screening of All Stoichiometric Inorganic Materials" *Chem* **1**, 617 (2016)](http://www.cell.com/chem/abstract/S2451-9294(16)30155-3)

[K. T. Butler et al, 
"Computational materials design of crystalline solids", *Chemical Society Reviews* (2016)](http://pubs.rsc.org/en/content/articlelanding/2016/cs/c5cs00841g)

[B. R. Pamplin, "A systematic method of deriving new semiconducting
compounds by structural analogy", *J. Phys. Chem. Solids*
**25**, 675 (1964)](http://www.sciencedirect.com/science/article/pii/0022369764901763)

[Mendeley "Materials Design" Reading List](https://www.mendeley.com/groups/8113991/materials-design/overview/)

Development notes
-----------------

Code style should comply with [PEP
8](http://www.python.org/dev/peps/pep-0008) where possible.
[Google's house
style](http://google-styleguide.googlecode.com/svn/trunk/pyguide.html)
is also helpful, including a good model for docstrings.
Please use comments liberally when adding nontrivial features, and
take the chance to clean up other people's code while looking at it.
 
The project was started Python 2.7.x, but has now been ported to Python 3. Please use new-style classes and string formatting.

Testing modules should be pass/fail and wrapped into **tests/test.py**.
Tests need to be run from the main directory (i.e. with `python tests/test.py`).
