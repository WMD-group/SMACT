smact.dopant_prediction module
==============================

The dopant prediction module facilitates high-throughput prediction of p-type and n-type dopants in multi-component solids. Candidates are filtered on accessible oxidation states, an n-type dopant being required to carry a higher charge than the host ion it replaces and a p-type dopant a lower one, with the sign of the site preserved. Survivors are scored using a data-mined substitution table or a species embedding, combined by default with a selectivity term. No ionic-radius criterion is applied; :mod:`smact.dopant_prediction.doper` explains why, and where to look if you do have a structure.

#.
   :class:`~smact.dopant_prediction.doper.Doper`, A class to identify possible dopants. :meth:`~smact.dopant_prediction.doper.Doper.get_dopants` returns a dictionary keyed by "n-type cation substitutions", "p-type cation substitutions", "n-type anion substitutions" and "p-type anion substitutions".


Submodules
----------

.. toctree::

    smact.dopant_prediction.doper
