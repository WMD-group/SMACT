Structure Prediction Module
===========================

The structure prediction module implements a minimalist
framework for high-throughput prediction of new compounds
based on species similarity indices. The typical workflow for
using this package is as follows:

#. Create a list of
   :class:`~smact.structure_prediction.structure.SmactStructure` s
   from a series of POSCAR files, the
   `Materials Project <https://www.materialsproject.org>`_
   or from a :class:`~smact.structure_prediction.database.StructureDB`.
#. Implement a :class:`~smact.structure_prediction.mutation.CationMutator`
   to perform species substitution probability calculations from a
   lambda table.
#. Feed the list of
   :class:`~smact.structure_prediction.structure.SmactStructure` s
   into the :class:`~smact.structure_prediction.mutation.CationMutator`
   and get a list of potential new
   :class:`~smact.structure_prediction.structure.SmactStructure` s.

Submodules
----------

.. toctree::

    smact.structure_prediction.database
    smact.structure_prediction.mutation
    smact.structure_prediction.structure
    smact.structure_prediction.probability_models
    smact.structure_prediction.utilities
