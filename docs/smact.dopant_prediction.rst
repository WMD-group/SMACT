Dopant Prediction Module
===========================

The dopant prediction module facilitates high-throughput prediction of p-type and n-type dopants in multi-component solids. The search and ranking process is based on electronic filters (e.g. accessible oxidation states) and chemical filters (e.g. difference in ionic radius).    

#. 
   :class:`~smact.dopant_prediction.doper`, A class to identify possible dopants. The output is the candidate species in the form of a dictionary with "n_type_cation", "p_type_cation", "n_type_anion", "p_type_anion", ranked by the corresponding elemental substitution probabilities.


Submodules
----------

.. toctree::

    smact.dopant_prediction.doper
 
