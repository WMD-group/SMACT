Structure Prediction Module
===========================

The structure prediction module implements a minimalist
framework for high-throughput prediction of new compounds
based on species similarity indices. The typical workflow for
using this package is as follows:

#. Create a database (:class:`~smact.structure_prediction.database.StructureDB`)
   of :class:`~smact.structure_prediction.structure.SmactStructure` s
   from the `Materials Project <https://www.materialsproject.org>`_.
#. Feed the database into a
   :class:`~smact.structure_prediction.prediction.StructurePredictor`,
   along with a list of chemical species,
   to predict potential likely structures
   by mutating structures in the database,
   such that they contain the given species.

Submodules
----------

.. toctree::

    smact.structure_prediction.prediction
    smact.structure_prediction.database
    smact.structure_prediction.mutation
    smact.structure_prediction.structure
    smact.structure_prediction.probability_models
    smact.structure_prediction.utilities
