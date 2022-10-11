.. _smact_module:

smact Python package
====================

The core module of :mod:`smact` contains classes which are used as
fundamental data types within the smact package, as well as several
utility functions.
Particular attention is drawn to :func:`smact.element_dictionary`,
which returns a dictionary of :class:`smact.Element` objects indexed
by their chemical symbols.
Generating this dictionary once and then performing lookups is
generally the fastest way of accessing element data while enumerating
possibilities.

.. automodule:: smact
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

.. toctree::

  smact.structure_prediction
  smact.dopant_prediction
  smact.properties
  smact.screening
  smact.oxidation_states
  smact.builder
  smact.distorter
  smact.lattice
  smact.lattice_parameters
  smact.data_loader
