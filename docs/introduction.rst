
Introduction
============

:mod:`smact` is a collection of tools and examples for "low-fi" screening of
potential semiconducting materials through the use of chemical
rules.

:mod:`smact` uses a combination of heuristics and models derived from data to
rapidly search large areas of chemical space. This combination of methods
allows :mod:`smact` to identify new materials for applications such as photovoltaics,
water splitting and thermoelectrics. 

Features of :mod:`smact` include:

- Chemical elements with associated properties
- Filters for oxidation states and charge balancing
- Structure predcition from chemical composition
- Composition probability prediction

Install
=======

The package is available *via* :code:`pip install smact`.

License and citation
====================

:mod:`smact` is distributed under an MIT license.

To cite the theory of :mod:`smact` please use:

- `Computational screening of all stoichiometric inorganic materials <https://www.sciencedirect.com/science/article/pii/S2451929416301553>`_

To cite the code of :mod:`smact` please use:

- `SMACT: Semiconducting Materials by Analogy and Chemical Theory <https://joss.theoj.org/papers/10.21105/joss.01361.pdf>`_

Studies using smact
===================

Read more about :mod:`smact` in our publications:

- `Computational screening of all stoichiometric inorganic materials <https://www.sciencedirect.com/science/article/pii/S2451929416301553>`_
- `Computer-aided design of metal chalcohalide semiconductors: from chemical composition to crystal structure <http://pubs.rsc.org/en/content/articlehtml/2017/sc/c7sc03961a>`_
- `Materials discovery by chemical analogy: role of oxidation states in structure prediction <http://pubs.rsc.org/en/content/articlehtml/2018/fd/c8fd00032h>`_

This approach is inspired by the work of Harrison [1]_ and
Pamplin [2]_. The work is an active project in the `Materials Design Group <http://wmd-group.github.io>`_.


We are also developing a set of Jupyter Notebook examples `here <https://github.com/WMD-group/SMACT/tree/master/examples>`_.

.. [1] http://www.worldcat.org/oclc/5170450 Harrison, W. A. *Electronic structure and the properties of solids: the physics of the chemical bond* (1980)

.. [2] http://dx.doi.org/10.1016/0022-3697(64)90176-3 Pamplin, B. R. *A systematic method of deriving new semiconducting compounds by structural analogy* J. Phys. Chem. Solids **7**, 675--684 (1964)
