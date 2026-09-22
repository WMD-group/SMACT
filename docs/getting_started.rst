
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

Alternatively, the latest version from the default branch of the Git repo can be installed using:

.. code::

    pip install git+https://github.com/WMD-group/SMACT.git

SMACT is also available from conda-forge:

.. code::

    conda install -c conda-forge smact

Optional functionality is provided by conda-forge feature packages. For example,
install the pre-trained property-prediction functionality with:

.. code::

    conda install -c conda-forge smact-property-prediction

The available feature packages are ``smact-mp``, ``smact-crystal-space``,
``smact-featurisers``, ``smact-visualisation``, ``smact-ml``, and
``smact-optional``. The latter matches the ``smact[optional]`` pip extra. To
install the dependencies required by all documented examples and tutorials:

.. code::

    conda install -c conda-forge smact-optional smact-property-prediction

These packages use conda-forge builds of compiled dependencies such as
``pytorch`` (the package is named ``torch`` on PyPI). Choose any GPU-specific
``pytorch`` configuration separately for your platform.

For developer installation, clone the repository and use `uv <https://docs.astral.sh/uv/>`_:

.. code::

    git clone https://github.com/wmd-group/smact.git
    cd smact
    uv sync --all-extras
    uv run pre-commit install

See ``CONTRIBUTING.md`` for full development workflow details.
