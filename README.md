SMACT
=====

Semiconducting Materials by Analogy and Computational Techniques

A collection of fast screening tools from elemental data.

Contents
--------

* **smact_core.py** 
* **smact_data.py** A collection of functions providing elemental properties.
* **smact_lattice.py** Given the sites, multiplicities and possible oxidation states at those sites, this reads from the database and generates all possible stoichiometeries.

Requirements
------------

The main language will be Python with Numpy, Scipy and Matplotlib.
The [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase) 
(ASE) is required for some components.

The [chemlab](http://chemlab.github.com/chemlab) project is not
currently used, but is considered ``friendly''; we will try to avoid
namespace clashes and it may be used for some features in the future.
Needed modifications will be made in [ajjackson's
fork](https://github.com/ajjackson/chemlab), 
but are expected to make it upstream fairly rapidly.



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
compounds by structural analogy", (1964) *J. Phys. Chem. Solids*
**25** pp. 675-684


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
