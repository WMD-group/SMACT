Dopant Prediction Module
===========================

The dopant prediction module facilitates high-throughput prediction of p-type and n-type dopants in multi-component solids. The search and ranking process is based on electronic filters (e.g. accessible oxidation states) and chemical filters (e.g. difference in ionic radius).    

#. 
   :class:`~smact.dopant_prediction.doper`, A class to identify possible n-type p-type dopants + give suggestions of original species given as a dictionary with "n_type_cation", "p_type_cation", "n_type_anion", "p_type_anion", and it calculates the corresponding substitution probabilities.


Submodules
----------

.. toctree::

    smact.dopant_prediction.doper
 
