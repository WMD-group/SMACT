SMACT
=====

Semiconducting Materials by Analogy and Computational Techniques

A collection of fast screening tools from elemental data.

Contents
--------

* **smact** library containing:
  * **core.py** 
  * **data.py** A collection of functions providing elemental properties.
  * **lattice.py** Given the sites, multiplicities and possible oxidation states at those sites, this reads from the database and generates all possible stoichiometeries.
  * **builder.py** Builds some common lattice structures, given the chemical composition.
  * **distorter.py** A collection of functions for enumerating and then substituting on inequivalent sites of a sub-lattice.

Requirements
------------

The main language is Python with Numpy, Scipy and Matplotlib.
The [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase) 
(ASE) is required for some components, as is [spglib](http://spglib.sourceforge.net).

The [chemlab](http://chemlab.github.com/chemlab) project is not
currently used, but is considered ``friendly''; we will try to avoid
namespace clashes and it may be used for some features in the future.
Needed modifications will be made in [ajjackson's
fork](https://github.com/ajjackson/chemlab), 
but are expected to make it upstream fairly rapidly.

Installation
------------

On a unix-like system, simply add the directory containing this README file
to your PYTHONPATH. e.g. in ~/.bashrc

    export PYTHONPATH="/home/username/src/smact:$PYTHONPATH"

License and attribution
-----------------------

SMACT is produced by the Walsh Materials Design group of the
Department of Chemistry at the University of Bath, UK.  Python code
and original data tables are licensed under the GNU General Public
License (GPL) v3.

The following files have their own licenses: **data/elements.txt ** is
from the OpenBabel project and licensed under the GPL v2, which is
included in the parent folder.

References
----------

B. R. Pamplin, "A systematic method of deriving new semiconducting
compounds by structural analogy", *J. Phys. Chem. Solids*
**25**, 675 (1964)

S. Chen, X. G. Gong, A. Walsh and S.-H. Wei, "Electronic structure and stability of quaternary chalcogenide semiconductors derived from cation cross-substitution of II-VI and I-III-VI2 compounds", *Physical Review B* **79**, 165211 (2009)

Development notes
-----------------

This project is currently in an early state, and is being developed by
members of the Walsh Materials Design group in the Department of
Chemistry at the University of Bath, UK.

Code style should comply with [PEP
8](http://www.python.org/dev/peps/pep-0008) where possible.
[Google's house
style](http://google-styleguide.googlecode.com/svn/trunk/pyguide.html)
is also helpful, including a good model for docstrings.
Please use comments liberally when adding nontrivial features, and
take the chance to clean up other people's code while looking at it.

Testing modules should be pass/fail and wrapped into **tests/test.py**.
Tests need to be run from the main directory (i.e. with `python tests/test.py`).


General To Do list:
-------------------
