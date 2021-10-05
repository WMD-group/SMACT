smact.properties module
=======================

A collection of tools for estimating useful properties.

The "electronegativity of a compound" computed with
:func:`compound_electroneg` is the rescaled geometric mean of
electronegativity used in Nethercot's recipe for estimating the
photoelectric threshold: [1]_

.. math:: \Phi^{AB} = 2.86(\chi_{A}\chi_{B})^{1/2} + E_{g} / 2.

In other words, the computed group
:math:`2.86(\chi_{A}\chi_{B})^{1/2}`
is the mid-gap energy and the VBM/CBM positions can be estimated by
subtracting/adding half of the band gap :math:`E_g`.
This is an extension Mulliken's electronegativity scale in which
:math:`\chi_{A} = (I_{A} + E_{A})/2` (where :math:`I` and :math:`E`
are respectively the ionisation potential and electron affinity.) [2]_

.. [1] Nethercot, A. H. (1974). *Phys. Rev. Lett.*, **33**, 1088–1091. http://dx.doi.org/10.1103/PhysRevLett.33.1088

.. [2] Mulliken, R. S. (1934). *J. Chem. Phys.*, **2**, 782. http://dx.doi.org/10.1063/1.1749394

.. automodule:: smact.properties
    :members:
    :undoc-members:
    :show-inheritance:
