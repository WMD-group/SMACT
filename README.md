SMACT
=====

Semiconducting Materials by Analogy and Computational Techniques

Contents
--------

* **smact_data.py** A collection of functions providing elemental properties. Draws on `chemlab`'s database (which in turn draws on openbabel)
* **smact_lattice.py** Given the sites, multiplicities and possible oxidation states at those sites, this reads from the database and generates all possible stoichiometeries.

Requirements
------------

The main language will be Python with Numpy, Scipy and Matplotlib.
The [chemlab](http://chemlab.github.com/chemlab) project is also used.
Modifications will be made in [ajjackson's fork](https://github.com/ajjackson/chemlab).

License and attribution
-----------------------

SMACT is produced by the Walsh Materials Design group of the Department of Chemistry at the University of Bath, UK. 
Python code and original data tables are licensed under the GNU General Public License (GPL) v3.

The following files have their own licenses:
**data/elements.txt ** is from the OpenBabel project and licensed under the GPL v2, which is included in the parent folder.

References
----------

B. R. Pamplin, "A systematic method of deriving new semiconducting compounds by structural analogy", (1964) *J. Phys. Chem. Solids* **25** pp. 675-684 
