
Getting Started
===============

============
Requirements
============

The main language is Python 3 and has been tested using Python 3.6+. Basic requirements are Numpy and Scipy.
The `Atomic Simulation Environment <https://wiki.fysik.dtu.dk/ase>`_
(ASE), `spglib <http://atztogo.github.io/spglib>`_,
and `pymatgen <http://pymatgen.org>`_ are also required for many components.

============
Installation
============

The latest stable release of SMACT can be installed via pip which will automatically setup other Python packages as required:

.. code::

    pip install smact

Alternatively, the latest master branch from the Git repo can be installed using:

.. code::

    pip install git+git://github.com/WMD-group/SMACT.git

Then ensure that the location of :mod:`smact` is on your PYTHONPATH.

For developer installation SMACT can be installed from a copy of the source repository (https://github.com/wmd-group/smact);
this will be preferred if using experimental code branches.

To clone the project from Github and make a local installation:

.. code::

    git clone https://github.com/wmd-group/smact.git
    cd smact
    pip install --user -e .

With -e pip will create links to the source folder so that that changes to the code will be immediately reflected on the PATH.
