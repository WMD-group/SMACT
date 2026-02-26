
Getting Started
===============

============
Requirements
============

The main language is Python 3 and has been tested using Python 3.11+.
Core dependencies include NumPy, SciPy, pandas, `pymatgen <http://pymatgen.org>`_,
`ASE <https://wiki.fysik.dtu.dk/ase>`_, and `spglib <http://atztogo.github.io/spglib>`_.
A full list is in ``pyproject.toml``.

============
Installation
============

The latest stable release of SMACT can be installed via pip, which will automatically setup other Python packages as required:

.. code::

    pip install smact

Alternatively, the latest master branch from the Git repo can be installed using:

.. code::

    pip install git+https://github.com/WMD-group/SMACT.git

For developer installation, clone the repository and use `uv <https://docs.astral.sh/uv/>`_:

.. code::

    git clone https://github.com/wmd-group/smact.git
    cd smact
    uv sync --all-extras
    uv run pre-commit install

See ``CONTRIBUTING.md`` for full development workflow details.
